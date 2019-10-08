plan = drake_plan(
  # annot = tibble(transcript = gencodev31$transcript_id, 
  #                gene_name = gencodev31$gene_name) %>%
  #   distinct() %>%
  #   filter(!is.na(transcript)),
  
  ## Read and process metadata:
  # found elsewhere that samples with a concentration below this level cluster together and differ substantially
  
  md =
    import_table(file=metadata_file,
                 bucket="memory_alpha",
                 FUN=read_excel, 
                 sheet = "main", 
                 col_types = c("numeric", 
                               rep("text",10),
                               "numeric",
                               rep("text",6),
                               rep("numeric", 3),
                               "text"),
                 trim_ws = TRUE,
                 na = "n/a",
                 .name_repair = ~ make_clean_names) %>%
    mutate(sample_name = make_clean_names(sample_name, case = "all_caps")) %>%
    filter(final_concentration_ng_ul > initial_concentration_threshold,
           project %in% (projects_to_include %||% unique(.data$project)),
           project %nin% projects_to_exclude,
           disease_class %in% (disease_classes_to_include %||% unique(.data$disease_class)),
           disease_class %nin% disease_classes_to_exclude) %>% 
    #Select the portions of the metadata that are useful:
    select(sample_name,
           study_group,
           disease_class,
           project,
           run_id,
           sample_alias,
           ord,
           sex,
           race_code,
           initial_concentration_ng_ul,
           final_concentration_ng_ul,
           rin) %>%
    mutate(project = factor(project),
           study_group = factor(study_group),
           disease_class = factor(disease_class),
           run_id = factor(run_id),
           initial_concentration_ng_ul = replace_na(initial_concentration_ng_ul, 200)),
  
  #Inspect the metadata:
  md_cat_data = inspect_cat(md),
  md_num_data = inspect_num(md),
  
  # Import count data
  #Raw FASTQs were processed using Salmon, which created pseudocounts.  
  #To import those, we need a named list (e.g. a dictionary) of the location
  #of all files ending in "h5" which each list member named for the sample it represents.
  
  tx_sample_names = dir(path = seq_file_directory,
                        pattern = "quant.sf.gz",
                        recursive = TRUE,
                        full.name = TRUE) %>% 
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>% 
    str_split(pattern = "/") %>% 
    map(function(x)`[[`(x,length(x)-1)) %>%
    # strip the sequencing run part ("_S156_L002") of the name
    str_split(pattern = '_') %>% 
    map_chr(`[[`,1) %>%
    make_clean_names(case = "all_caps"),
  
  tx_files = dir(path = seq_file_directory,
                 pattern = "quant.sf.gz",
                 recursive = TRUE,
                 full.name = TRUE) %>% 
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
    `names<-`(tx_sample_names) %>%
    `[`(!is.na(match(names(.), md$sample_name))),
  
  final_md = md[which(md[['sample_name']] %in% names(tx_files)),] %>%
    column_to_rownames('sample_name'),
  samples = tx_files[rownames(final_md)],
  counts = tximport(samples,
                    type = "salmon",
                    txIn = TRUE,
                    txOut = FALSE,
                    tx2gene = annot,
                    importer = fread),
  
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normCts = counts$counts/counts$length,
  
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib = calcNormFactors(normCts) * colSums(normCts),
  
  normMat = (counts$length/exp(rowMeans(log(counts$length)))) %>%
    sweep(x = .,
          MARGIN = 2,
          STATS = eff.lib,
          FUN = "*") %>%
    log(),
  
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  
  pre_qc_dge =
    DGEList(counts=counts$counts,
            samples = final_md,
            group = final_md[,comparison_grouping_variable],
            lib.size = eff.lib,
            remove.zeros = TRUE) %>%
    scaleOffset(y = .,
                offset = normMat[rownames(.$counts),]) %>%
    `[`(filterByExpr(y = .,
                     group = final_md[,comparison_grouping_variable],
                     keep.lib.sizes=FALSE,
                     min.count = 0,
                     min.total.count = 10),),
  
  preliminary_design = model.matrix(study_design,
                         final_md[colnames(pre_qc_dge),
                                  c("initial_concentration_ng_ul",
                                    "run_id",
                                    "disease_class")]) %>%
    `colnames<-`(c("Intercept",
                   colnames(.)[2:ncol(.)])),
  
  outlier_qc = remove_outliers(object = pre_qc_dge,
                                zscore_cutoff = pc1_zscore_threshold,
                                design = preliminary_design),
  
  pca_qc = outlier_qc$pca,
  removed_outliers = outlier_qc$removed,
  
  design = model.matrix(study_design,
                         final_md[colnames(outlier_qc$count_object),
                                  c("initial_concentration_ng_ul",
                                    "run_id",
                                    "disease_class")]) %>%
    `colnames<-`(c("Intercept",
                   colnames(.)[2:ncol(.)])),
  
  # contr.matrix = makeContrasts(experimental_group = colnames(design)[ncol(design)],
  #                               levels=colnames(design)),
  
  # Creating a DGEList object for use in edgeR.
  dge_qc = outlier_qc$count_object %>% 
    scaleOffset(y = ., 
                offset = normMat[rownames(.$counts),colnames(.$counts)]) %>%
    `[`(filterByExpr(y = .,
                     group = final_md[,comparison_grouping_variable],
                     keep.lib.sizes=FALSE,
                     min.count = 0,
                     min.total.count = 10), 
    ) %>%
    calcNormFactors(object = .) %>%
    estimateDisp(y = ., design = design, robust = TRUE),
  
  fit = glmQLFit(y = dge_qc,
                 design = dge_qc$design),
  qlf = glmQLFTest(fit, coef=ncol(dge_qc$design)),
  
  res = topTags(qlf, n = Inf) %>%
    `$`('table') %>%
    as_tibble(rownames = "gene") %>% 
    arrange(desc(logFC)),
  
  v = voom(counts = dge_qc, design = dge_qc$design),
  
  # Minus the noise, actual differences and similarities are now apparent.
  sample_dists = v$E %>% t() %>% dist(),
  sampleDistMatrix = as.matrix(sample_dists),
  sample_dendrogram = sample_dists %>% hclust() %>% as.dendrogram(),
  
  
  annotation_info = as.data.frame(dge_qc$samples)[,c("disease_class",
                                                             "study_group",
                                                             "run_id",
                                                             "sex",
                                                             "initial_concentration_ng_ul",
                                                             "final_concentration_ng_ul",
                                                             "rin")],
  
  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_res = prcomp_irlba(x = v$E) %>%
    `[[`("rotation") %>%
    as_tibble() %>%
    mutate(sample = colnames(v$E)) %>% 
    inner_join(as_tibble(v$targets, 
                         rownames="sample")),
  
  #fig.width=12, fig.height=6
  pca_plot = pca_res %>%
    ggplot(aes(x = PC1,
               y = PC2,
               color = disease_class)) +
    geom_point(alpha = 0.3) + 
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot(),
  
  #fig.width=12, fig.height=6
  pca_plot2 = pca_res %>%
    ggplot(aes(x = PC1,
               y = PC2,
               color = run_id)) +
    geom_point(alpha = 0.3) + 
    #geom_mark_hull(aes(fill=run_id)) +
    scale_color_paletteer_d(package="ggsci",
                            palette = "category20_d3") +
    scale_fill_paletteer_d(package="ggsci",
                           palette = "category20_d3") +
    labs(color = "run_id") +
    theme_cowplot(),
  
  umap_results = umap(t(v$E),
                      n_threads = detectCores(),
                      n_sgd_threads = detectCores(),
                      verbose = TRUE) %>% 
    as_tibble(.name_repair = "unique") %>%
    `colnames<-`(c("umap_1",
                  "umap_2")) %>%
    mutate(sample = colnames(v$E)) %>% 
    inner_join(as_tibble(v$targets,
                         rownames="sample")),
  
  umap_plot = umap_results %>% 
    ggplot(aes(x = umap_1,
               y = umap_2,
               color = disease_class)) +
    geom_point(alpha=0.5) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(color = "disease_class") +
    theme_cowplot(),
  
  umap_plot2 = umap_results %>% 
    ggplot(aes(x = umap_1, 
               y = umap_2,
               color = run_id)) +
    geom_point(alpha=0.5) +
    scale_color_paletteer_d(package="ggsci",
                            palette = "category20_d3") +
    scale_fill_paletteer_d(package="ggsci",
                           palette = "category20_d3") +
    theme_cowplot(),
  
  down_table =
    res %>%
    filter(FDR <= 0.05) %>%
    mutate(logFC = -logFC) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(25, logFC) %>%
    arrange(desc(logFC)),
  
  up_table = 
    res %>%
    filter(FDR <= 0.05) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(25, logFC) %>%
    arrange(desc(logFC)),
  
  degs = 
    res %>% 
    filter(FDR < 0.05),
  
  top_up = degs %>%
    filter(logFC > 0) %>%
    top_n(25, logFC) %>%
    pull(gene),
  
  top_down = degs %>%
    filter(logFC < 0) %>%
    top_n(25, -logFC) %>%
    pull(gene),
  
  ISGs = intersect(c("STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
                     "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1", 
                     "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3", "IFITM1", "IRF9", "ISG15", 
                     "IFI6", "IFIT3", "USP18", "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", 
                     "IRF3", "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27", "GBP2", 
                     "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2", "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"),
                   rownames(v$E)),
  
  module_tbl = module_list %>%
    enframe() %>%
    unnest(cols = "value") %>%
    rename(module = "name",
           gene = "value") %>%
    inner_join(as_tibble(module_annot,
                         rownames="module")),
  
  module_ISGs = module_tbl %>%
    dplyr::group_by(gene) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(type == "Interferon",
                  gene %in% rownames(v$E)) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene"),
  
  sle_module_scores = scoreEigengenes(object = v$E,
                                      module_list = module_list,
                                      score_func = 'rsvd') %>%
    mutate(sample = colnames(dge_qc)) %>%
    select(sample, everything()),
  
  ldg_module_scores = scoreEigengenes(object = v$E,
                                      module_list = ldg_modules,
                                      score_func = 'rsvd') %>%
    mutate(sample = colnames(dge_qc)) %>%
    select(sample, everything()),
  
  annotated_modules = module_annot %>% as_tibble(rownames="module") %>% filter(type != "Undetermined"),
  
  type_pal = paletteer_d(package="ggsci",
                         palette = "category20_d3",
                         n = length(unique(annotated_modules$type))) %>%
    `names<-`(unique(annotated_modules$type)),
  
  # we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
  run_groups =
    dge_qc$samples %>%
    as_tibble() %>%
    pull(run_id) %>%
    unique(),
 
  run_id_color_set =
    colorRampPalette(
      brewer.pal(12, "Paired"))(length(run_groups)) %>%
    `names<-`(run_groups),
  
  chr_pal = c("Y" = "#0000FF", "X" = "#FF0000"),
  
  sex_pal = c("Male" = "#0000FF", "Female" = "#FF0000"),
  
  comparison_grouping_variable_colors = c("#000000",
                                          "#FF9999") %>% `names<-`(c(control_group, experimental_group)),
  
  group_pal =
    list(
      comparison_grouping_variable_colors,
      run_id_color_set,
      type_pal,
      chr_pal,
      sex_pal
      ) %>% `names<-`(c(comparison_grouping_variable, "run_id", "type", "chr", "sex")),
  
  module_scores =
    sle_module_scores %>%
    inner_join(ldg_module_scores),
  
  annotated_module_scores =
    module_scores %>%
    select(sample, 
           one_of(annotated_modules$module),
           ldg_a,
           ldg_b) %>%
    column_to_rownames("sample"),

  viral_transcripts =
    annot %>%
    filter(!str_detect(string = transcript,
                       pattern = "^ENST")) %>%
    pull(gene_name) %>%
    intersect(rownames(v$E)),
  
  viral_exprs =
    v$E[viral_transcripts,] %>%
    t() %>%
    as_tibble(rownames="sample"),
  
  ifn_modules =
    annotated_modules %>%
    filter(type == "Interferon") %>%
    pull(module),
  
  inflame_modules =
    annotated_modules %>%
    filter(type == "Inflammation") %>%
    pull(module),
  
  ifn_scores =
    module_scores %>%
    select(sample,
           one_of(ifn_modules)),
  
  inflammation_scores =
    module_scores %>%
    select(sample,
           one_of(inflame_modules)),
  
  ldg_scores =
    module_scores %>%
    select(sample,
           one_of(names(ldg_modules))),
  
  ifn_scores_with_viral =
    inner_join(viral_exprs, ifn_scores),
  
  inflammation_scores_with_viral =
    inner_join(viral_exprs, inflammation_scores),
    
  report = rmarkdown::render(
    knitr_in("report.rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE)
)
