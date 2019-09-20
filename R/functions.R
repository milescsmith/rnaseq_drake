`%nin%` <- compose(`!`, `%in%`)

deduplicate_samples <- function(md, samples){
  if (nrow(get_dupes(md, sample_name)) > 0){
    deduplicated_md = md %>% 
      filter(sample_name %in% names(samples)) %>%
      mutate(sample_name = make_clean_names(string = sample_name, 
                                            case = "all_caps"))
    deduplicated_samples = `names<-`(samples, make_clean_names(string = names(samples),
                                                                 case = "all_caps"))  
  } else {
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples))
    deduplicated_samples = samples
  }
  
  return(list(md = deduplicated_md,
           samples = deduplicated_samples))
}

remove_outliers <- function(dds, zscore_cutoff){
  dds  <- estimateSizeFactors(dds,
                              locfun = shorth,
                              type = "poscounts")
  vsd <- assay(vst(dds))
  pca_res = prcomp_irlba(vsd)[['rotation']] %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd),
           zscore = abs((PC1 - mean(PC1))/sd(PC1))) %>%
    inner_join(as_tibble(colData(dds), rownames = "sample"))
  
  outliers <- pca_res %>% filter(zscore >= zscore_cutoff) %>% pull(sample)
  if (length(outliers > 0)){
    dds <- dds[,colnames(dds) %nin% outliers] 
  }
  return(list(dds = dds, pca = pca_res, removed = outliers))
}


top_DEGs <- function(results_obj, number_DEGs = 100, padj_cutoff = 0.05){
  
  results_tbl = 
    results_obj %>% 
    as_tibble(rownames = "gene") %>%
    filter(padj < padj_cutoff)
  
  top_up = results_tbl %>%
    filter(log2FoldChange > 0) %>%
    top_n(number_DEGs, log2FoldChange) %>%
    pull(gene)
  
  top_down = results_tbl %>%
    filter(log2FoldChange < 0) %>%
    top_n(number_DEGs, -log2FoldChange) %>%
    pull(gene)
  
  return(list(up = top_up, 
              down = top_down))
}

DEG_tables <- function(results_obj, number_DEGs = 25, padj_cutoff = 0.05){
  
  down_table = 
    results_obj %>%
    as_tibble(rownames = "gene") %>%
    filter(padj <= padj_cutoff) %>%
    mutate(log2FoldChange = -log2FoldChange) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(number_DEGs, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
  
  up_table = 
    results_obj %>%
    as_tibble(rownames = "gene") %>%
    filter(padj <= padj_cutoff) %>%
    mutate_at(vars(-gene), list(~signif(., 2))) %>%
    top_n(number_DEGs, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
  
  return(list(up = up_table,
              down = down_table))
}

split_process <- function(mother_obj, selected_celltype, selected_ethnicity, zscore_cutoff, comparison_var){
  
    comparison_var <- enexpr(comparison_var)
    
    child_obj_pre_qc = colData(mother_obj) %>% 
      as_tibble(rownames = "sample") %>%
      filter(celltype == selected_celltype,
             ethnicity == selected_ethnicity) %>%
      pull(sample) %>%
      `[`(x = mother_obj, j = .)
    
    filtered_samples = remove_outliers(dds = child_obj_pre_qc,
                                       zscore_cutoff = zscore_cutoff)
    pca_qc = filtered_samples$pca
    design(filtered_samples$dds) <- as.formula(paste("~", quo_text(comparison_var)))
    dds = DESeq(filtered_samples$dds, parallel = TRUE)
    
    compare_combos = combn(factor(mixedsort(unique(dds[[quo_text(comparison_var)]]))), 2)
    res_combinations = map(1:ncol(compare_combos), function(x){
      results(
        object = dds,
        contrast = c(quo_text(comparison_var),
                     compare_combos[2,x],
                     compare_combos[1,x]),
        alpha = 0.05,
        parallel = TRUE)
    })
    
    names(res_combinations) <- map_chr(1:ncol(compare_combos), 
                                       function(x){str_glue("{compare_combos[1,x]}_vs_{compare_combos[2,x]}")})
    
    DEG_tables_res = map(res_combinations, DEG_tables)
    top_DEG_res = map(res_combinations, top_DEGs)
    
    return(list(dds = dds,
                pca_qc = pca_qc,
                res = res_combinations,
                tables = DEG_tables_res,
                top_degs = top_DEG_res,
                vsd = assay(vst(dds))))
}
