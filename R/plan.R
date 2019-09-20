plan = drake_plan(
  # annot = tibble(transcript = gencodev31$transcript_id, 
  #                gene_name = gencodev31$gene_name) %>%
  #   distinct() %>%
  #   filter(!is.na(transcript)),
  
  ## Read and process metadata:
  # found elsewhere that samples with a concentration below this level cluster together and differ substantially
  
  md =
    read_csv(file = "/home/milo/datasets/bulk_preclinical/original_patient_md.csv") %>%
    clean_names() %>%
    mutate(ethnicity = factor(ethnicity),
           disease_class = factor(class),
           celltype = factor(celltype),
           sample = make_clean_names(sample, case = "all_caps")) %>%
    select(-class, -rin, -conc, -cell_count),
  
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
    # When str_split splits a string, it makes everything before the matching pattern into an element of the returned list
    # even if there is nothing before the split - you just get an empty element
    # thus, the seventh element matches '012210101_S156_L002'
    map_chr(function(x)`[[`(x,length(x)-1)) %>%
    make_clean_names(case = "all_caps"),
  
  tx_files = dir(path = seq_file_directory,
                 pattern = "quant.sf.gz",
                 recursive = TRUE,
                 full.name = TRUE) %>% 
    grep(pattern = "Undetermined|NONE", invert = TRUE, value = TRUE) %>%
    `names<-`(tx_sample_names) %>%
    `[`(!is.na(match(names(.), md$sample))),
  
  #Most data structures that support row or column names cannot tolerate duplicates.  Are there any duplicate samples that need to be fixed?
  # md_dupes = md %>% filter(sample_name %in% tx_sample_names) %>% get_dupes(sample_name),
  # 
  # dedupe = deduplicate_samples(md, tx_files),
  # deduplicated_md = dedupe$md,
  # deduplicated_tx_files = dedupe$sample,
  
  # Read in count files
  # final_md = deduplicated_md[which(deduplicated_md[['sample_name']] %in% names(deduplicated_tx_files)),] %>%
  #   column_to_rownames('sample_name'),
  final_md = md[which(md[['sample']] %in% names(tx_files)),] %>%
    column_to_rownames('sample'),
  samples = tx_files[rownames(final_md)],
  counts = tximport(samples,
                    type = "salmon",
                    txIn = TRUE,
                    txOut = FALSE,
                    tx2gene = annot),
  
  # Process the data using DESeq2
  # The `DESeq` function run several subfunctions that take 
  # care of normalization, dispersion calculations, model fitting,
  # and differential expression analysis.  This can take quite some 
  # time, especially as the study design grows more complex.
  dds_import = DESeqDataSetFromTximport(txi = counts,
                                        colData = final_md,
                                        design = study_design), 
  dds_filtered = dds_import %>%
    `[`(rowSums(counts(.)) > 1, ) %>%
    `[`(grep(pattern = "^RNA5",
             x = rownames(.),
             invert = TRUE,
             value = TRUE),),
  
  ## Sample QC filtering
  # Remove samples that have a PC1 Z-score > 3. This matches what I was doing visually, but is vastly quicker.
  outlier_qc = remove_outliers(dds = dds_filtered,
                               zscore_cutoff = 3),
  dds_qc = outlier_qc$dds,
  pca_qc = outlier_qc$pca,
  removed_outliers = outlier_qc$removed,
  
  # I would just use the DESeq() function, but running each contitutent separately makes it easier to recover from
  # failure or to observe progress
  
  dds_processed = DESeq(dds_qc,
                        parallel = TRUE),
  vsd = vst(dds_processed),
  vsd_exprs = assay(vsd),
  
  sample_dists = vsd_exprs %>% t() %>% dist(),
  #fig.width=12, fig.height=9
  sampleDistMatrix = as.matrix(sample_dists),
  
  annotation_info = as.data.frame(colData(dds_processed))[,c("disease_class",
                                                             "ethnicity",
                                                             "celltype")],
  
  # Minus the noise, actual differences and similarities are now apparent.
  sample_dendrogram = sample_dists %>% hclust() %>% as.dendrogram(),
  
  # Variation can also be examined in reduced dimensional space by PCA or UMAP:
  pca_res = prcomp_irlba(x = vsd_exprs) %>%
    `[[`("rotation") %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd_exprs)) %>% 
    inner_join(as_tibble(colData(dds_processed), 
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
  umap_results = umap(t(vsd_exprs),
                      n_threads = detectCores(),
                      n_sgd_threads = detectCores(),
                      verbose = TRUE) %>% 
    as_tibble(.name_repair = "unique") %>%
    `colnames<-`(c("umap_1",
                   "umap_2")) %>%
    mutate(sample = colnames(vsd_exprs)) %>% 
    inner_join(as_tibble(colData(dds_processed),
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
  
  
  ea_t = split_process(mother_obj = dds_processed,
                       selected_celltype = "T",
                       selected_ethnicity = "EA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  ea_b = split_process(mother_obj = dds_processed,
                       selected_celltype = "B",
                       selected_ethnicity = "EA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  ea_m = split_process(mother_obj = dds_processed,
                       selected_celltype = "M",
                       selected_ethnicity = "EA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  aa_t = split_process(mother_obj = dds_processed,
                       selected_celltype = "T",
                       selected_ethnicity = "AA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  aa_b = split_process(mother_obj = dds_processed,
                       selected_celltype = "B",
                       selected_ethnicity = "AA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  aa_m = split_process(mother_obj = dds_processed,
                       selected_celltype = "M",
                       selected_ethnicity = "AA",
                       zscore_cutoff = 3,
                       comparison_var = disease_class),
  
  #fig.width=10, fig.height=9
  ISGs = intersect(c("STAT1", "ADAR", "ABCE1", "RNASEL", "TYK2", "IFNAR1",
                     "IFNB1", "STAT2", "IFNAR2", "JAK1", "SAMHD1", 
                     "SOCS1", "SOCS3", "STAT1", "ISG20", "IFITM3", "IFITM1", "IRF9", "ISG15", 
                     "IFI6", "IFIT3", "USP18", "IP6K2", "PSMB8", "IFIT1", "IRF4", "IRF5", "IRF1", 
                     "IRF3", "IRF6", "IRF8", "IRF2", "IRF7", "IFITM2", "XAF1", "IFI27", "GBP2", 
                     "RSAD2", "MX2", "MX1", "IFIT2", "IFI35", "BST2", "OAS1", "OASL", "OAS2", "OAS3", "PTPN6", "PTPN11"),
                   rownames(vsd_exprs)),
  
  module_tbl = module_list %>%
    enframe() %>%
    unnest(cols = "value") %>%
    rename(name = "module",
           value = "gene") %>%
    inner_join(as_tibble(module_annot,
                         rownames="module")),
  
  module_ISGs = module_tbl %>%
    dplyr::group_by(gene) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::filter(type == "Interferon",
                  gene %in% rownames(vsd_exprs)) %>%
    dplyr::select(module, gene) %>%
    arrange(module) %>%
    column_to_rownames("gene"),
  
  dds_with_scores = scoreEigengenes(dds_processed, module_list = module_list),
  
  annotated_modules = module_annot %>% as_tibble(rownames="module") %>% filter(type != "Undetermined"),
  
  type_pal = paletteer_d(package="ggsci",
                         palette = "category20_d3",
                         n = length(unique(annotated_modules$type))) %>%
    `names<-`(unique(annotated_modules$type)),
  
  # we manually setup the palettes for pheatmap because letting it automatically pick the colors results in terrible choices
  
  class_pal = paletteer_d("ggthemes", "calc")[1:3] %>% `names<-`(c("NEG","POS","SLE")),
  ethnicity_pal = paletteer_d("ggthemes", "calc")[c(4,11)] %>% `names<-`(c("AA", "EA")),
  celltype_pal = paletteer_d("ggthemes", "calc")[8:10] %>% `names<-`(c("T","B","M")),
  pathway_pal = paletteer_d("awtools", "ppalette")[1:7] %>% `names<-`(unique(sgl$pathway)),
  group_pal =
    list(
      class_pal,
      ethnicity_pal,
      celltype_pal,
      type_pal,
      pathway_pal
    ) %>% `names<-`(c("disease_class", "ethnicity", "celltype", "type", "pathway")),
  
  module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample") %>%
    select(sample, 
           matches("^M[[:digit:]]+\\.")) %>%
    column_to_rownames("sample"),
  
  annotated_module_scores =
    colData(dds_with_scores) %>%
    as_tibble(rownames="sample") %>%
    select(sample, 
           one_of(annotated_modules$module)) %>%
    column_to_rownames("sample"),
  
  ifn_modules =
    annotated_modules %>%
    filter(type == "Interferon") %>%
    pull(module),
  
  inflame_modules =
    annotated_modules %>%
    filter(type == "Inflammation") %>%
    pull(module),
  
  ifn_scores =
    colData(dds_with_scores)[,ifn_modules] %>%
    as.data.frame() %>%
    as_tibble(rownames="sample"),
  
  inflammation_scores =
    colData(dds_with_scores)[,inflame_modules] %>%
    as.data.frame() %>%
    as_tibble(rownames="sample"),
  
  report = rmarkdown::render(
    knitr_in("report.rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE)
)
