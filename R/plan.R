md =
  read_excel(
    path = metadata_file,
    sheet = "main",
    trim_ws = TRUE,
    na = "n/a"
    ) %>%
  clean_names() %>%
  rename(
    sample_name = nova_seq_sample_id,
    ethnicity = race_code) %>%
  select(
    -ord,
    -pos,
    -i7_index,
    -index,
    -i5_index,
    -index2,
    -correction_needed,
    -sample_alias,
    -mess,
    -abc_or_mess_or_control
    ) %>%
  mutate(
    sample_name = janitor::make_clean_names(sample_name, case = "all_caps"),
    disease_class =
      str_remove(
        string =
          recode(
            disease_class,
            "Unaffected Control" = "control"
            ),
        pattern = " Patient"
      ),
    ethnicity =
      str_remove_all(
        string = ethnicity,
        pattern = "[\\[\\]]"
      ),
    initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200),
    across(
      .cols = c(project, run_id, study, sex, ethnicity, disease_class),
      .fns = as_factor
      ),
    across(
      .cols = c(age, contains("concentration")),
      .fns = as.numeric
      )
  ) %>%
  filter(
    initial_concentration_ng_ul > initial_concentration_threshold,
    project %in% (projects_to_include %||% unique(.data$project)),
    project %nin% projects_to_exclude,
    disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
    disease_class %nin% disease_classes_to_exclude
    ) %>%
  mutate(
    across(
      .cols = where(is.factor),
      .fns = ~ fct_drop(.x)
    )
  ) %>%
  distinct()

tx_files =
  find_sequencing_samples(
    seq_file_directory = seq_file_directory,
    source_names = md[["sample_name"]],
    alignment_source = "salmon"
    )

final_md =
  md %>%
  filter(
    sample_name %in% names(tx_files),
    str_detect(
      string = sample_name,
      pattern = "_2$",
      negate = TRUE
    )
  ) %>%
  mutate(
    disease_class = fct_relevel(disease_class, {{control_group}}),
    run_id = fct_drop(run_id)
  ) %>%
  column_to_rownames('sample_name')

#Inspect the metadata:
md_cat_data = inspectdf::inspect_cat(final_md)

md_num_data = inspectdf::inspect_num(final_md)

pruned_samples = tx_files[rownames(final_md)]

experimental_results =
  process_counts(
    count_files                  = tx_files,
    sample_metadata              = final_md,
    study_design                 = study_design,
    comparison_grouping_variable = comparison_grouping_variable,
    aligner                      = "salmon",
    minimum_gene_count           = 10,
    pc1_zscore_threshold         = pc1_zscore_threshold,
    pc2_zscore_threshold         = pc2_zscore_threshold,
    BPPARAM                      = BPPARAM,
    sva_num                      = num_sva,
    use_combat                   = FALSE,
    method                       = analysis_library
  )

  transformed_counts = experimental_results[["transformed_counts"]]

  banchereau_module_scores =
    tirosh_score_modules(
      expr_obj    = transformed_counts,
      module_list = banchereau_modules
      ) %>%
    as_tibble(rownames = "sample_name")

  ldg_module_scores =
    tirosh_score_modules(
      expr_obj    = transformed_counts,
      module_list = ldg_modules
      ) %>%
    as_tibble(rownames = "sample_name")

  metasig_scores =
    tirosh_score_modules(
      expr_obj    = transformed_counts,
      module_list = metasignature_module
      ) %>%
    as_tibble(rownames = "sample_name")

  module_scores =
    inner_join(
      banchereau_module_scores,
      ldg_module_scores
      ) %>%
    inner_join(metasig_scores) %>%
    column_to_rownames(var = "sample_name")

  counts_obj_with_scores =
    add_sample_metadata(
      object = experimental_results[["dataset"]],
      md = module_scores
      )

  #
  # scoreEigengenes(object = dds_processed,
  #                 module_list = module_list,
  #                 score_func = 'rsvd') %>%
  # scoreEigengenes(object = .,
  #                 module_list = ldg_modules,
  #                 score_func = 'rsvd') %>%
  # scoreEigengenes(object = .,
  #                 module_list = metasignature_module,
  #                 score_func = 'rsvd'),

  sample_dists =
    transformed_counts %>%
    t() %>%
    parallelDist::parallelDist()

  #fig.width=12, fig.height=9
  sampleDistMatrix = as.matrix(sample_dists)

  sample_dendrogram =
    sample_dists %>%
    hclust() %>%
    as.dendrogram()

  annotated_modules =
    module_annotation %>%
    filter(type != "Undetermined")

  annotated_mod_list =
    annotated_modules %>%
    mutate(
      module_type =
        paste(
          module,
          type,
          sep = " - "
          )
    ) %>%
    select(-type) %>%
    deframe()

  annotated_module_scores =
    module_scores %>%
    select(one_of(annotated_modules$module))

  sample_cluster_info =
    para_ident_clusters(
      annotated_module_scores,
      K.max = 20
      )

  clusters = sample_cluster_info$clusters %>% mutate(cluster = as_factor(cluster))

  counts_obj_with_scores_clusters =
    add_sample_metadata(
      object = counts_obj_with_scores,
      md = clusters
      )

  annotation_info =
    retrieve_metadata(
      object = counts_obj_with_scores_clusters,
      columns = c("disease_class", "sex","cluster","sample_name")
      ) %>%
    column_to_rownames(var = "sample_name")

  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_results = irlba::prcomp_irlba(x = transformed_counts) %>%
    pluck("rotation") %>%
    as_tibble() %>%
    mutate(sample_name = colnames(transformed_counts)) %>%
    inner_join(retrieve_metadata(object = counts_obj_with_scores_clusters))

  umap_results =
    umap(
      t(transformed_counts),
      n_threads = detectCores(),
      n_sgd_threads = detectCores(),
      verbose = TRUE,
      n_components = 3
    ) %>%
    as_tibble(.name_repair = "unique") %>%
    set_names(
      c(
        "umap_1",
        "umap_2",
        "umap_3"
        )
      ) %>%
    mutate(sample_name = colnames(transformed_counts)) %>%
    inner_join(retrieve_metadata(object = counts_obj_with_scores_clusters))

  comparison_results_list = experimental_results[["comparisons"]]

  res = experimental_results[["degs"]]

  down_tables = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        !is.na(padj) & padj <= 0.05,
        log2FoldChange < 0
      ) %>%
      mutate(
        log2FoldChange = -log2FoldChange
      ) %>%
      mutate_at(
        .vars = vars(-gene),
        .funs = list(~signif(x = ., digits =  2))
      ) %>%
      top_n(
        n = 25,
        wt = log2FoldChange
      ) %>%
      arrange(
        desc(
          log2FoldChange
        )
      )
  }) %>%
    set_names(
      nm = map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    ) %>%
    keep(~ nrow(.x) > 0)

  up_tables = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        !is.na(padj) & padj <= 0.05,
        log2FoldChange > 0
      ) %>%
      mutate_at(
        vars(-gene),
        list(~signif(., 2)
        )
      ) %>%
      top_n(
        n = 25,
        wt = log2FoldChange
      ) %>%
      arrange(
        desc(
          log2FoldChange
        )
      )
  }) %>%
    set_names(
      nm = map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    ) %>%
    keep(~ nrow(.x) > 0)

  degs = map(seq_along(res), function(i){
    res[[i]] %>%
      as_tibble(
        rownames = "gene"
      ) %>%
      filter(
        padj < 0.05
      ) %>%
      filter(
        abs(
          log2FoldChange
        ) >= 0.5
      ) %>%
      pull(gene)
  }) %>%
    set_names(
      map_chr(
        .x = comparison_results_list,
        .f = str_remove,
        pattern = str_glue("{comparison_grouping_variable}_")
      )
    )

  deg_class = enframe(degs) %>%
    unnest(cols = c(value)) %>%
    group_by(value) %>%
    mutate(count = n()) %>%
    mutate(name = case_when(count == 1 ~ name,
                            count > 1 ~ "multiple")) %>%
    select(-count) %>%
    distinct() %>%
    column_to_rownames(var = "value") %>%
    set_names("comparison")

  deg_means =
    transformed_counts %>%
    t() %>%
    as_tibble(rownames = "name") %>%
    select(name, one_of(rownames(deg_class))) %>%
    pivot_longer(-name,
                 names_to = "gene",
                 values_to = "expr") %>%
    left_join(
      as_tibble(
        colData(dds_processed),
        rownames = "name"
        ) %>%
        select(
          name,
          sex,
          ethnicity,
          disease_class
          )
      ) %>%
    group_by(gene,
             disease_class) %>%
    summarise(avg = mean(expr)) %>%
    pivot_wider(names_from = gene,
                values_from = avg) %>%
    column_to_rownames("disease_class")

  #fig.width=10, fig.height=9
  ISGs = intersect(c("STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
                     "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1",
                     "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3",
                     "IFITM1", "IRF9", "ISG15", "IFI6", "IFIT3", "USP18",
                     "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", "IRF3",
                     "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27",
                     "GBP2", "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2",
                     "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"),
                   rownames(transformed_counts))

  module_tbl =
    c(
      banchereau_modules,
      ldg_modules,
      metasignature_module
      ) %>%
    enframe() %>%
    unnest(cols = "value") %>%
    rename(module = name,
           gene = value) %>%
    inner_join(module_annotation)

  module_ISGs = module_tbl %>%
    dplyr::group_by(gene) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(type == "Interferon",
                  gene %in% rownames(transformed_counts)) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene")

  #### PALETTES ####
  # we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
  type_pal =
    paletteer_d(
      "ggsci::category20_d3",
      n = length(levels(annotated_modules$type))
      ) %>%
    as.character() %>%
    set_names(levels(annotated_modules$type)) %>%
    magrittr::inset2("Undetermined", "#000000")

  chr_pal = c("Y" = "#E41A1C",
              "X" = "#377EB8")

  sex_pal = c("Male" = "coral3",
              "Female" = "azure2",
              "unk" = "#333333")

  cluster_pal =
    ifelse(
      test = length(levels(clusters$cluster)) > 12,
      yes = list(
        colorRampPalette(
          paletteer_d(
            palette = "ggthemes::calc",
            n = 12
            )
          )(
            length(
              levels(
                clusters$cluster
                )
              )
            )
        ),
      no = list(
        paletteer_d(
          palette = "ggthemes::calc",
          n = length(
            levels(
              clusters$cluster
              )
            )
          )
        )
      ) %>%
    unlist() %>%
    as.character() %>%
    set_names(levels(clusters$cluster))

  project_pal =
    colorRampPalette(
      brewer.pal(9, "Set1"))(length(levels(annotation_info$project))) %>%
    set_names(levels(annotation_info$project))

  number_disease_classes =
    length(
      unique(
        annotation_info$disease_class
        )
      )

  disease_class_pal =
    if_else(
      number_disease_classes > 2,
      list(RColorBrewer::brewer.pal(number_disease_classes, "Set1")),
      list(c("black", "grey75"))
    ) %>%
    unlist() %>%
    set_names(
      unique(
        annotation_info$disease_class
        )
      )

  # cell_type_pal =
  #   c("#ffa600", "#0062cc", "#008a71") %>%
  #   set_names(levels(annotation_info$cell_type))

  comparison_pal =
    oaColors::oaPalette(
      length(
        unique(deg_class$comparison)
        )
      ) %>%
    set_names(unique(deg_class$comparison))

  group_pal =
    list(
      type_pal,
      chr_pal,
      sex_pal,
      cluster_pal,
      project_pal,
      disease_class_pal,
      comparison_pal ) %>% #,
      # cell_type_pal) %>%
    set_names(c(
      "type",
      "chr",
      "sex",
      "cluster",
      "project",
      "disease_class",
      "comparison")) #,
      #"cell_type"))

  top_vars =
    matrixStats::rowSds(transformed_counts) %>%
    set_names(rownames(transformed_counts)) %>%
    enframe() %>%
    top_n(20000, value) %>%
    pull(name)

  vsd_top = transformed_counts[top_vars,] %>%
    t()

  sft = pickSoftThreshold(
    data = vsd_top,
    powerVector =
      c(
        seq(10),
        seq(
          from = 12,
          to = 20,
          by = 2)
        ),
    verbose = 5)

  #### WGCNA ####
  wgcna_modules =
    blockwiseModules(
      datExpr = vsd_top,
      power = sft$powerEstimate,
      maxBlockSize = 20000,
      mergeCutHeight = 0.2,
      minModuleSize = 20,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 3,
      detectCutHeight = 0.995,
      TOMDenom = "min",
      networkType = "signed hybrid",
      reassignThreshold = 1e-6
    )

  wgcna_module_genes =
    wgcna_modules$colors %>%
    enframe(
      name = "gene",
      value = "module"
      )

  wgcna_hub_genes =
    chooseTopHubInEachModule(
      datExpr = vsd_top,
      colorh = wgcna_modules$colors,
      power = 4,
      type = "signed hybrid"
      )

  wgcna_scores = wgcna_modules$MEs %>%
    as_tibble(rownames="sample_name") %>%
    select(-MEgrey) %>%
    left_join(as_tibble(annotation_info,
                        rownames="sample_name"))

  wgcna_cluster_split =
    initial_split(
      data = wgcna_scores %>%
        mutate(disease_class = fct_drop(disease_class)) %>%
        select(
          cluster,
          starts_with("ME")
          ),
      prop = 0.75,
      strata = "cluster"
      )

  wgcna_cluster_train = training(wgcna_cluster_split)

  wgcna_cluster_test = testing(wgcna_cluster_split)

  wgcna_cluster_rf_cv =
    train(
      cluster ~ .,
      method = "parRF",
      data = wgcna_cluster_train,
      trControl =
        trainControl(
          method = "repeatedcv",
          number = 10,
          repeats = 10,
          search = "grid",
          allowParallel = TRUE
        ),
      importance=T
    )

  wgcna_cluster_rf_cv_varImp =
    varImp(
      object = wgcna_cluster_rf_cv,
      scale = FALSE,
      importance = TRUE
      )

  wgcna_disease_class_split =
    initial_split(
      data = wgcna_scores %>%
        # filter(disease_class != "LP") %>%
        mutate(disease_class = fct_drop(disease_class)) %>%
        select(
          disease_class,
          starts_with("ME")
          ),
      prop = 0.75,
      strata = "disease_class"
      )

  wgcna_disease_class_train = training(wgcna_disease_class_split)

  wgcna_disease_class_test = testing(wgcna_disease_class_split)

  wgcna_disease_class_rf_cv =
    train(
      form = disease_class ~ .,
      method = "parRF",
      data = wgcna_disease_class_train,
      trControl =
        trainControl(
          method = "repeatedcv",
          number = 20,
          repeats = 20,
          search = "grid",
          allowParallel = TRUE
          )
      )

  wgcna_disease_class_rf_cv_varImp =
    varImp(
      object = wgcna_disease_class_rf_cv,
      scale = FALSE,
      importance=TRUE
      )

  #### module_classification ####
  module_scores_with_md =
    module_scores %>%
    as_tibble(rownames="sample_name") %>%
    left_join(as_tibble(annotation_info,
                        rownames="sample_name"))

  module_cluster_split =
    initial_split(
      data = module_scores_with_md %>%
        mutate(disease_class = fct_drop(disease_class)) %>%
        select(
          cluster,
          one_of(names(banchereau_modules))
          ),
      prop = 0.75,
      strata = "cluster"
      )

  module_cluster_train = training(module_cluster_split)

  module_cluster_test = testing(module_cluster_split)

  module_cluster_rf_cv =
    train(
      form = cluster ~ .,
      method = "parRF",
      data = module_cluster_train,
      trControl =
        trainControl(
          method = "repeatedcv",
          number = 10,
          repeats = 10,
          search = "grid",
          allowParallel = TRUE
        ),
      importance = TRUE
    )

  module_cluster_rf_cv_varImp =
    varImp(
      object = module_cluster_rf_cv,
      scale = FALSE,
      importance = TRUE
    )

  module_disease_class_split =
    initial_split(
      data = module_scores_with_md %>%
        mutate(disease_class = fct_drop(disease_class)) %>%
        select(
          disease_class,
          one_of(names(banchereau_modules))
        ),
      prop = 0.75,
      strata = "disease_class"
    )

  module_disease_class_train = training(module_disease_class_split)

  module_disease_class_test = testing(module_disease_class_split)

  module_disease_class_rf_cv =
    train(
      form = disease_class ~ .,
      method = "parRF",
      data = module_disease_class_train,
      trControl =
        trainControl(
          method = "repeatedcv",
          number = 10,
          repeats = 10,
          search = "grid",
          allowParallel = TRUE
        ),
      importance = TRUE
    )

  module_disease_class_rf_cv_varImp =
    varImp(
      object = module_disease_class_rf_cv,
      scale = FALSE,
      importance = TRUE
    )

  #### WGCNA-based GSEA ####
  filtered_wgcna_module_genes =
    wgcna_module_genes %>%
    mutate(hugo = checkGeneSymbols(gene)[["Suggested.Symbol"]]) %>%
    filter(!is.na(hugo))

  MEenriched =
    future_map_dfr(
      .x = unique(filtered_wgcna_module_genes$module),
      .f = function(i){
        filtered_wgcna_module_genes %>%
          filter(module == i) %>%
          pull(hugo) %>%
          enricher(gene = .,
                   TERM2GENE = c5)%>%
          {if(!is.null(.)) mutate(slot(.,'result'), module = i) } #if_else does not work because it tries to evaulate both the true and false results, which doesn't work with `slot(NULL, "results")`
      })

  MEplotting = MEenriched %>%
    filter(p.adjust < 0.05,
           module != "grey") %>%
    mutate(
      GeneRatio = map_dbl(
        .x = GeneRatio,
        .f = function(i){
          j = str_split(i, "/") %>%
            magrittr::extract2(1) %>%
            as.double()
          j[[1]]/j[[2]]
          }),
      ID = str_replace_all(string = ID, pattern = "_", replacement = " "),
      module = paste0("ME", module)
      ) %>%
    group_by(module) %>%
    top_n(5, GeneRatio) %>%
    sample_n(size = 5,
             replace = TRUE) %>%
    distinct() %>%
    ungroup() %>%
    arrange(module, GeneRatio) %>%
    mutate(order = row_number())

  viral_transcripts =
    annot %>%
    filter(
      !str_detect(
        string = transcript,
        pattern = "^ENST"
        )
      ) %>%
    pull(gene_name) %>%
    intersect(rownames(transformed_counts))

  detected_viral_transcripts =
    counts(counts_obj_with_scores) %>%
    t() %>%
    as_tibble(rownames = "sample_name") %>%
    select(
      sample_name,
      one_of(viral_transcripts)
      ) %>%
    pivot_longer(
      -sample_name,
      names_to = "transcript",
      values_to = "counts"
      ) %>%
    group_by(transcript) %>%
    summarise(total_counts = sum(counts)) %>%
    filter(total_counts > 0) %>%
    pull(transcript)

  viral_exprs =
    transformed_counts[detected_viral_transcripts,] %>%
    t() %>%
    as_tibble(rownames = "sample_name")

  ifn_modules =
    annotated_modules %>%
    filter(type == "Interferon") %>%
    pull(module)

  inflame_modules =
    annotated_modules %>%
    filter(type == "Inflammation") %>%
    pull(module)

  module_scores_with_viral =
    retrieve_metadata(counts_obj_with_scores_clusters) %>%
    select(
      sample_name,
      disease_class,
      # cell_type,
      one_of(names(banchereau_modules)),
      one_of(names(ldg_modules)),
      names(metasignature_module)
    ) %>%
    inner_join(wgcna_scores) %>%
    inner_join(viral_exprs) %>%
    inner_join(clusters)

  annotated_module_scores_with_cluster_class =
    module_scores_with_viral %>%
    mutate(
      cluster = as_factor(cluster)
      # cell_type = as_factor(cell_type)
      ) %>%
    select(
      cluster,
      # cell_type,
      disease_class,
      one_of(annotated_modules$module)
    )

  renamed_annotated_module_scores =
    names(annotated_module_scores_with_cluster_class) %>%
    drake_recode(thing_to_unquote_splice = annotated_mod_list) %>%
    set_names(
      nm = .,
      x = annotated_module_scores_with_cluster_class
    )

  annotated_module_scores_pivot =
    renamed_annotated_module_scores %>%
    pivot_longer(
      cols = starts_with("M"),
      names_to = "module",
      values_to = "score"
    )

  annotated_module_stats_by_cluster =
    annotated_module_scores_pivot %>%
    group_by(module) %>%
    wilcox_test(
      score ~ cluster,
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = annotated_module_scores_pivot,
      group_var = module,
      compare_value = score
    )

  # annotated_module_stats_by_cell_type =
  #   annotated_module_scores_pivot %>%
  #   group_by(module) %>%
  #   wilcox_test(
  #     score ~ cell_type,
  #     p.adjust.method = "BH"
  #   ) %>%
  #   grouped_add_xy_positions(
  #     stats_tbl = .,
  #     data_tbl = annotated_module_scores_pivot,
  #     group_var = module,
  #     compare_value = score
  #   )

  annotated_module_stats_by_disease =
    annotated_module_scores_pivot %>%
    group_by(module) %>%
    wilcox_test(
      score ~ disease_class,
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = annotated_module_scores_pivot,
      group_var = module,
      compare_value = score
    )

  module_scores_pivot =
    module_scores_with_viral %>%
    mutate(
      cluster = as_factor(cluster),
      disease_class = as_factor(disease_class)
    ) %>%
    select(
      sample_name,
      cluster,
      disease_class,
      matches("^M[[:digit:]]+"),
      mg,
      starts_with("ldg")
    ) %>%
    pivot_longer(
      cols = c(
        matches("^M[[:digit:]]+"),
        mg,
        starts_with("ldg")
      ),
      names_to = "module",
      values_to = "score"
    )

  module_stats_by_cluster =
    module_scores_pivot %>%
    group_by(module) %>%
    wilcox_test(
      score ~ cluster,
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = module_scores_pivot,
      group_var = module,
      compare_value = score
    )

  module_stats_by_disease =
    module_scores_pivot %>%
    group_by(module) %>%
    wilcox_test(
      score ~ disease_class,
      p.adjust.method = "BH"
    ) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = module_scores_pivot,
      group_var = module,
      compare_value = score
    )

  module_scores_with_viral_by_cluster =
    module_scores_with_viral %>%
    select(cluster,
           matches("^ME")) %>%
    pivot_longer(-cluster,
                 names_to = "module",
                 values_to = "score") %>%
    mutate(cluster = as_factor(cluster)
    )

  module_scores_with_viral_by_cluster_stats =
    module_scores_with_viral_by_cluster %>%
    group_by(module) %>%
    wilcox_test(score ~ cluster) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = module_scores_with_viral_by_cluster,
      group_var = module,
      compare_value = score
    )

  # module_scores_with_viral_by_cell_type =
  #   module_scores_with_viral %>%
  #   select(cell_type,
  #          matches("^ME")) %>%
  #   pivot_longer(-cell_type,
  #                names_to = "module",
  #                values_to = "score") %>%
  #   mutate(cell_type = as_factor(cell_type)
  #   ),

  # module_scores_with_viral_by_cell_type_stats =
  #   module_scores_with_viral_by_cell_type %>%
  #   group_by(module) %>%
  #   wilcox_test(score ~ cell_type) %>%
  #   grouped_add_xy_positions(
  #     stats_tbl = .,
  #     data_tbl = module_scores_with_viral_by_cell_type,
  #     group_var = module,
  #     compare_value = score
  #   ),

  module_scores_with_viral_by_disease =
    module_scores_with_viral %>%
    select(disease_class,
           matches("^ME")) %>%
    pivot_longer(-disease_class,
                 names_to = "module",
                 values_to = "score") %>%
    mutate(disease_class = as_factor(disease_class)
    )

  module_scores_with_viral_by_disease_stats =
    module_scores_with_viral_by_disease %>%
    group_by(module) %>%
    wilcox_test(score ~ disease_class) %>%
    grouped_add_xy_positions(
      stats_tbl = .,
      data_tbl = module_scores_with_viral_by_disease,
      group_var = module,
      compare_value = score
    )

  transformed_counts %>%
    as.data.frame() %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/filtered_normalized_stabilized_expression.csv")
      )

  retrieve_metadata(counts_obj_with_scores_clusters) %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/filtered_metadata.csv")
      )

  counts(counts_obj_with_scores) %>%
    as.data.frame() %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/filtered_transcript_counts.csv")
      )

  res %>%
    as.data.frame() %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/log_fold_changes.csv")
      )

  wgcna_modules$MEs %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/wgcna_eigengene_scores.csv")
      )

  save_wgcna_module_genes =
    wgcna_module_genes %>%
    write_csv(path = file_out("results/wgcna_genes.csv"))

  retrieve_metadata(counts_obj_with_scores_clusters) %>%
    select(
      sample_name,
      matches("^M[[:digit:]]+"),
      one_of(names(ldg_modules))
      ) %>%
    as.data.frame() %>%
    data.table::fwrite(
      row.names = TRUE,
      file = file_out("results/module_scores.csv")
      )

  disp_plot = plot_dispersion_estimate(counts_obj_with_scores)

  module_scores_for_correlation =
    module_scores_with_viral %>%
    select(
      sample_name, `M4.10`, `M4.11`, `M4.1`, `M4.15`, `M3.6`, `M8.46`, `M4.14`,
      `M5.15`, `ldg_1.1`, `ldg_2.1`, `M8.83`, `M3.2`, `M4.2`, `M4.6`, `M4.13`,
      `M5.1`, `M5.7`, `M7.1`, `M1.2`, `M3.4`, `M5.12`, `M6.6`, `M2.2`, `M3.3`,
      `M3.5`, `M4.7`, `M6.11`,`M6.16`, `M9.42`, `M6.13`, `M5.6`, `M5.10`,
      `M6.2`, `M6.12`, `M4.3`, `M4.5`, `M5.9`, `M1.1`, `M2.3`,`M3.1`, `M6.18`
    ) %>%
    set_names(
      nm = c(
        "sample_name", "M4.10_Bcell", "M4.11_PlasmaCells", "M4.1_Tcells",
        "M4.15_Tcells", "M3.6_Cytotoxic_NK", "M8.46_Cytotoxic_NK",
        "M4.14_Monocytes", "M5.15_Neutrophils","LDG1.1", "LDG2.1",
        "M8.83_ImmuneResponses", "M3.2_Inflammation", "M4.2_Inflammation",
        "M4.6_Inflammation", "M4.13_Inflammation", "M5.1_Inflammation",
        "M5.7_Inflammation","M7.1_Inflammation", "M1.2_Interferon",
        "M3.4_Interferon", "M5.12_Interferon","M6.6_Apoptosis_Survival",
        "M2.2_CellCycle", "M3.3_CellCycle", "M3.5_CellCycle",
        "M4.7_CellCycle", "M6.11_CellCycle", "M6.16_CellCycle",
        "M9.42_CellCycle", "M6.13_CellDeath", "M5.6_Mitochondrial",
        "M5.10_Mitochondrial", "M6.2_Mitochondrial","M6.12_Mitochondrial",
        "M4.3_ProteinSynthesis", "M4.5_ProteinSynthesis",
        "M5.9_ProteinSynthesis", "M1.1_Platelets", "M2.3_Erythrocytes",
        "M3.1_Erythrocytes","M6.18_Erythrocytes"
      )
    )

  module_score_pivot =
    module_scores_for_correlation %>%
    pivot_longer(
      -sample_name,
      names_to = "module",
      values_to = "score"
    )

  module_score_correlation =
    module_scores_for_correlation %>%
    column_to_rownames("sample_name") %>%
    cor_mat(method = "spearman")

  report =
    rmarkdown::render(
      input = knitr_in("markdown/report.rmd"),
      output_file = file_out("results/report.html"),
      output_dir = "results",
      quiet = TRUE,
      params = "ask"
      )

  # report_pdf = rmarkdown::render(
  #   input = knitr_in("markdown/report_pdf.rmd"),
  #   output_file = file_out("results/report.pdf"),
  #   output_dir = "results",
  #   output_format = "all",
  #   quiet = TRUE),

  qc_report =
    rmarkdown::render(
      input = knitr_in("markdown/qc_report.rmd"),
      output_file = file_out("results/qc_report.html"),
      output_dir = "results",
      quiet = TRUE,
      params = "ask"
      )

  # supplemental_report =
  #   rmarkdown::render(
  #     input = knitr_in("markdown/supplemental_report.rmd"),
  #     output_file = file_out("results/supplemental_report.html"),
  #     output_dir = "results",
  #     output_format = "all",
  #     quiet = TRUE
  #     ),

  # manuscript_figures = rmarkdown::render(
  #   input = knitr_in("markdown/for_manuscript.rmd"),
  #   output_file = file_out("results/for_manuscript.pdf"),
  #   output_dir = "results",
  #   quiet = TRUE),

  # save(
  #   annotation_info,
  #   final_md,
  #   group_pal,
  #   module_scores_with_viral,
  #   ISGs,
  #   pca_results,
  #   degs,
  #   umap_results,
  #   transformed_counts,
  #   wgcna_modules,
  #   file = file_out("results/for_shiny_vis.RData")
  #   )
