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
