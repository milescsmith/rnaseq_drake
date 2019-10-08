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


#' Remove outliers
#' 
#' Perform PCA on a dataset and return one in which samples with a PC1
#' zscore greater than a given cutoff are removed
#'
#' @param object 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
remove_outliers <- function(object, ...){
  UseMethod("remove_outliers")
}

#' @rdname remove_outliers
#' @method remove_outliers DGEList
#' @importFrom irlba prcomp_irlba
#' @importFrom limma voom
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return list
remove_outliers.DGEList <- function(object, zscore_cutoff, design){
  v <- voom(object, design)
  pca_res = prcomp_irlba(v$E)[['rotation']] %>%
    as_tibble() %>%
    mutate(sample = colnames(v$E),
           zscore = abs((PC1 - mean(PC1))/sd(PC1))) %>%
    inner_join(as_tibble(v$targets, rownames = "sample"))
  
  outliers <- pca_res %>% filter(zscore >= zscore_cutoff) %>% pull(sample)
  if (length(outliers > 0)){
    object <- object[,colnames(object) %nin% outliers] 
  }
  return(list(count_object = object, pca = pca_res, removed = outliers))
}


#' @rdname remove_outliers
#' @method remove_outliers DESeqDataSet
#' @importFrom irlba prcomp_irlba
#' @importFrom DESeq2 vst estimateSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return DESeqDataSet
remove_outliers.DESeqDataSet <- function(object, zscore_cutoff){
  object  <- estimateSizeFactors(object,
                              locfun = shorth,
                              type = "poscounts")
  vsd <- assay(vst(object))
  pca_res = prcomp_irlba(vsd)[['rotation']] %>%
    as_tibble() %>%
    mutate(sample = colnames(vsd),
           zscore = abs((PC1 - mean(PC1))/sd(PC1))) %>%
    inner_join(as_tibble(colData(object), rownames = "sample"))
  
  outliers <- pca_res %>% filter(zscore >= zscore_cutoff) %>% pull(sample)
  if (length(outliers > 0)){
    object <- object[,colnames(object) %nin% outliers] 
  }
  return(list(count_object = object, pca = pca_res, removed = outliers))
}
