`%nin%` <- purrr::compose(`!`, `%in%`)
select <- dplyr::select

deduplicate_samples <- function(md, samples){
  if (nrow(get_dupes(md, sample_name)) > 0){
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples)) %>%
      mutate(sample_name = make_clean_names(string = sample_name,
                                            case = "all_caps"))
    deduplicated_samples = set_names(
      x = samples,
      nm = make_clean_names(
        string = names(samples),
        case = "all_caps"
      )
    )
  } else {
    deduplicated_md = md %>%
      filter(sample_name %in% names(samples))
    deduplicated_samples = samples
  }

  list(
    md = deduplicated_md,
    samples = deduplicated_samples
  )
}


#' @title process_counts
#'
#' @description Read in RNAseq counts and perform QC,
#' normalization, modeling, and DEG analysis
#'
#' @param ... Parameters to pass along to the actual functions
#' that will process the RNAseq data
#' @param method Process the data using limma, edgeR, or DESeq2?
#'
#' @return
#' @export
#'
#' @examples
process_counts <- function(..., method) {
  if (method == "limma"){
    process_counts.limma(...)
  } else if (method == "DESeq2") {
    process_counts.deseq2(...)
  } else if (method == "edgeR") {
    process_counts.edgeR(...)
  }
}

process_counts.limma <-
  function(
    count_files,
    sample_metadata,
    study_design,
    # batch_variable,
    comparison_grouping_variable,
    aligner                      = "salmon",
    minimum_gene_count           = 1,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold.        = 2,
    BPPARAM.                     = BPPARAM,
    sva_num                      = 2,
    use_combat                   = FALSE
  ){
    common_names <-
      intersect(
        x = names(count_files),
        y = rownames(sample_metadata)
      )

    count_files <- magrittr::extract(count_files, common_names)

    counts <-
      tximport(
        files    = count_files,
        type     = "salmon",
        txIn     = TRUE,
        txOut    = FALSE,
        tx2gene  = annot,
        importer = fread
      )

    message("Correcting for effective library sizes")
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normCts <- counts[["counts"]]/counts[["length"]]

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)

    # Multiply each gene by the library size for each sample and take the log
    # (sweep applies a function either rowwise or column wise; STATS is a vector
    # equal in length to the chosen axis to use as an argument to the applied function)
    normMat <-
      sweep(x      = counts[["length"]]/exp(rowMeans(log(counts[["length"]]))),
            MARGIN = 2,
            STATS  = eff.lib,
            FUN    = "*"
      ) %>%
      log()

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.

    pre_qc_dge <-
      DGEList(
        counts       = magrittr::extract(counts[["counts"]], ,rownames(sample_metadata)),
        samples      = sample_metadata,
        group        = sample_metadata[[comparison_grouping_variable]],
        lib.size     = eff.lib,
        remove.zeros = TRUE
      ) %>%
      scaleOffset(
        y = .,
        offset = normMat[rownames(.[["counts"]]),]
      ) %>%
      magrittr::extract(
        filterByExpr(
          y               = .,
          group           = sample_metadata[[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    preliminary_design <-
      model.matrix(
        object = study_design,
        data   = magrittr::extract(
          sample_metadata,
          colnames(pre_qc_dge),
          c("final_concentration_ng_ul",
            batch_variable,
            comparison_grouping_variable)
        )
      ) %>%
      magrittr::set_colnames(
        c("Intercept", colnames(.)[2:ncol(.)])
      )

    message("Removing outliers")
    outlier_qc <-
      remove_outliers(
        object            = pre_qc_dge,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold,
        design            = preliminary_design
      )

    design =
      model.matrix(
        object = study_design,
        data   =
          magrittr::extract(
            sample_metadata,
            colnames(outlier_qc$count_object),
            c("final_concentration_ng_ul",
              batch_variable,
              comparison_grouping_variable)
          )
      ) %>%
      set_colnames(
        c("Intercept",
          colnames(.)[2:ncol(.)])
      )


    # Creating a DGEList object for use in edgeR.
    pre_sva_dge <-
      scaleOffset(
        y                = outlier_qc[["count_object"]],
        offset           =
          normMat[
            rownames(outlier_qc[["count_object"]][["counts"]]),
            colnames(outlier_qc[["count_object"]][["counts"]])
          ]
      )

    pre_sva_dge <-
      magrittr::extract(
        pre_sva_dge,
        filterByExpr(
          y               = pre_sva_dge,
          group           = sample_metadata[[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    sva_res <-
      calc_sva(
        object       = pre_sva_dge,
        model_design = study_design,
        n.sva        = num_sva
      )

    post_qc_dge <-
      calcNormFactors(sva_res[["dge"]])

    groups <- sva_res[["dge"]][["samples"]][[comparison_grouping_variable]]

    mm <- model.matrix(~0 + groups)

    voom_exprs <-
      voomWithQualityWeights(
        counts    = post_qc_dge,
        design    = mm,
        plot      = TRUE,
        save.plot = TRUE
      )

    fit <-
      lmFit(
        object = voom_exprs,
        design = mm,
        method = "robust"
      )

    comparisons <-
      as_tibble(
        combinations(
          n = length(unique(groups)),
          r = 2,
          v = unique(groups)
        )
      ) %>%
      set_names(
        nm = c("V1","V2")
      ) %>%
      transmute(
        name = paste0(V2, "_vs_", V1),
        compare = paste(V2, "-", V1)
      ) %>%
      filter(
        str_detect(
          string = compare,
          pattern = control_group
        )
      ) %>%
      deframe()

    contr_matrix <-
      makeContrasts(
        contrasts = comparisons,
        levels = unique(groups)
      )

    colnames(fit$coefficients) <-
      str_remove(
        string  = colnames(fit$coefficients),
        pattern = "groups"
      )

    treat_res <-
      contrasts.fit(
        fit = fit,
        contrasts =  contr_matrix
      ) %>%
      treat()

    res = map(colnames(contr_matrix), function(i) {

      topTreat(fit = treat_res, coef = i) %>%
        as_tibble(rownames = "gene") %>%
        arrange(desc(logFC)) %>%
        rename(
          baseMean       = AveExpr,
          log2FoldChange = logFC,
          pvalue         = P.Value,
          padj           = adj.P.Val
        )
    }) %>%
      set_names(
        nm = names(comparisons)
      )

    list(
      raw_counts         = counts[["counts"]],
      normalized_counts  = cpm(post_qc_dge, log = TRUE),
      transformed_counts = voom_exprs[["E"]],
      outlier_samples.   = outlier_qc[["removed"]],
      qc_pca             = outlier_qc[["pca"]],
      degs               = res,
      dataset            = post_qc_dge,
      comparisons        = comparisons
    )
  }

process_counts.edgeR <-
  function(
    count_files,
    sample_metadata,
    study_design,
    # batch_variable,
    comparison_grouping_variable,
    aligner                      = "salmon",
    minimum_gene_count           = 1,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold.        = 2,
    BPPARAM.                     = BPPARAM,
    sva_num                      = 2,
    use_combat                   = FALSE
  ){
    common_names <-
      intersect(
        x = names(count_files),
        y = rownames(sample_metadata)
        )

    count_files <- magrittr::extract(count_files, common_names)

    counts <-
      tximport(
        files    = count_files,
        type     = "salmon",
        txIn     = TRUE,
        txOut    = FALSE,
        tx2gene  = annot,
        importer = vroom
      )

    message("Correcting for effective library sizes")
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normCts <- counts[["counts"]]/counts[["length"]]

    # Computing effective library sizes from scaled counts, to account for
    # composition biases between samples.
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)

    # Multiply each gene by the library size for each sample and take the log
    # (sweep applies a function either rowwise or column wise; STATS is a vector
    # equal in length to the chosen axis to use as an argument to the applied function)
    normMat <-
      sweep(x      = counts[["length"]]/exp(rowMeans(log(counts[["length"]]))),
            MARGIN = 2,
            STATS  = eff.lib,
            FUN    = "*"
      ) %>%
      log()

    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.

    pre_qc_dge <-
      DGEList(
        counts       = magrittr::extract(counts[["counts"]], ,rownames(sample_metadata)),
        samples      = sample_metadata,
        group        = sample_metadata[[comparison_grouping_variable]],
        lib.size     = eff.lib,
        remove.zeros = TRUE
      ) %>%
      scaleOffset(
        y = .,
        offset = normMat[rownames(.[["counts"]]),]
      ) %>%
      magrittr::extract(
        filterByExpr(
          y               = .,
          group           = sample_metadata[[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    preliminary_design <-
      model.matrix(
        object = study_design,
        sample_metadata[colnames(pre_qc_dge),
                        c("final_concentration_ng_ul",
                          "run_id",
                          "disease_class")
        ]
      ) %>%
      magrittr::set_colnames(
        c("Intercept", colnames(.)[2:ncol(.)])
      )

    message("Removing outliers")
    outlier_qc <-
      remove_outliers(
        object            = pre_qc_dge,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold,
        design            = preliminary_design
      )

    design =
      model.matrix(
        object = study_design,
        sample_metadata[colnames(outlier_qc$count_object),
                        c("initial_concentration_ng_ul",
                          batch_variable,
                          "disease_class")
        ]
      ) %>%
      set_colnames(
        c("Intercept",
          colnames(.)[2:ncol(.)])
      )


    # Creating a DGEList object for use in edgeR.

    pre_sva_dge <-
      scaleOffset(
        y                = outlier_qc[["count_object"]],
        offset           =
          normMat[
            rownames(outlier_qc[["count_object"]][["counts"]]),
            colnames(outlier_qc[["count_object"]][["counts"]])
          ]
      )

    pre_sva_dge <-
      magrittr::extract(
        pre_sva_dge,
        filterByExpr(
          y               = pre_sva_dge,
          group           = sample_metadata[[comparison_grouping_variable]],
          keep.lib.sizes  = FALSE,
          min.count       = minimum_gene_count,
          min.total.count = 10
        ),
      )

    sva_res <-
      calc_sva(
        object       = pre_sva_dge,
        model_design = study_design,
        n.sva        = num_sva
      )

    post_qc_dge <-
      calcNormFactors(
        object = sva_res[["dge"]]
      )

    post_qc_design <-
      model.matrix(
        object = sva_res[["design"]],
        data   = post_qc_dge[["samples"]]
      )

    post_qc_dge <-
      estimateDisp(
        y      = post_qc_dge,
        design = post_qc_design,
        robust = TRUE
      )

    fit =
      glmQLFit(
        y      = post_qc_dge,
        design = post_qc_design
      )

    comparison_results_list <-
      colnames(post_qc_design) %>%
      keep(
        str_detect(
          string  = .,
          pattern = comparison_grouping_variable
        )
      )

    res = map(comparison_results_list, function(i) {
      qlf =
        glmQLFTest(
          glmfit = fit,
          coef   = i
        )

      res = topTags(qlf, n = Inf) %>%
        magrittr::use_series('table') %>%
        as_tibble(rownames = "gene") %>%
        arrange(desc(logFC)) %>%
        rename(
          baseMean       = logCPM,
          log2FoldChange = logFC,
          pvalue         = PValue,
          padj           = FDR
        )
    }) %>%
      set_names(
        nm = map_chr(
          .x      = comparison_results_list,
          .f      = str_remove,
          pattern = comparison_grouping_variable
        ) %>%
          paste0("_vs_", control_group)
      )

    v <-
      voomWithQualityWeights(
        counts = post_qc_dge,
        design = post_qc_design
      )

    list(
      raw_counts         = counts(counts[["counts"]]),
      normalized_counts  = lcpm(post_qc_dge, log = TRUE),
      transformed_counts = v[["E"]],
      outlier_samples.   = outlier_qc[["removed"]],
      qc_pca             = outlier_qc[["pca"]],
      degs               = res,
      dataset            = post_qc_dge,
      comparisons        = comparison_results_list
    )
  }

process_counts.deseq2 <-
  function(
    count_files,
    sample_metadata,
    study_design,
    batch_variable,
    comparison_grouping_variable,
    aligner                      = "salmon",
    minimum_gene_count           = 1,
    pc1_zscore_threshold         = 2,
    pc2_zscore_threshold         = 2,
    BPPARAM                      = BPPARAM,
    use_combat                   = FALSE
  ){
    counts <-
      tximport(
        files    = count_files,
        type     = aligner,
        txIn     = TRUE,
        txOut    = FALSE,
        tx2gene  = annot,
        importer = fread
      )

    dds_import <-
      DESeqDataSetFromTximport(
        txi     = counts,
        colData = sample_metadata,
        design  = study_design
      )

    if (isTRUE(use_combat)){
      corrected_counts <-
        ComBat_seq(
          counts = counts(dds_import),
          batch  = fct_drop(colData(dds_import)[[batch_variable]]),
          group  = colData(dds_import)[[comparison_grouping_variable]]
        ) %>%
        `storage.mode<-`("integer")

      dds_import <-
        DESeqDataSetFromMatrix(
          countData = corrected_counts,
          colData   = colData(dds_import),
          design    = study_design
        )
    }

    dds_filtered <-
      dds_import %>%
      magrittr::extract(rowSums(counts(.)) > 1, ) %>%
      magrittr::extract(
        grep(
          pattern = "^RNA5",
          x       = rownames(.),
          invert  = TRUE,
          value   = TRUE
        )
      )

    ## Sample QC filtering
    # Remove samples that have a PC1 Z-score > 3. This matches what I was doing visually, but is vastly quicker.
    outlier_qc <-
      remove_outliers(
        count_obj         = dds_filtered,
        pc1_zscore_cutoff = pc1_zscore_threshold,
        pc2_zscore_cutoff = pc2_zscore_threshold
      )

    dds <-
      DESeq(
        object   = outlier_qc$dds,
        parallel = TRUE,
        BPPARAM  = BPPARAM
      )

    sva_res <- calc_sva(dds = dds, model_design = "disease_class", n.sva = num_sva)
    sva_graph_data <- sva_res$sva

    vsd <- vst(sva_res$dds)

    comparison_results_list <-
      resultsNames(object = sva_res$dds) %>%
      keep(
        str_detect(
          string  = .,
          pattern = comparison_grouping_variable
        )
      )

    res <-
      map(
        .x = comparison_results_list,
        .f = function(i) {
          lfcShrink(
            dds      = dds_processed,
            coef     = i,
            parallel = TRUE,
            type     = "apeglm")
        }
      ) %>%
      set_names(
        map_chr(
          .x      = comparison_results_list,
          .f      = str_remove,
          pattern = paste0(comparison_grouping_variable, "_")
        )
      )


    list(
      raw_counts         = counts(sva_res$dds),
      normalized_counts  = counts(sva_res$dds, normalized = TRUE),
      transformed_counts = assay(vsd),
      outlier_samples    = outlier_qc$removed,
      qc_pca             = outlier_qc$pca,
      sva_graph_data     = sva_res$sva,
      degs               = res,
      dataset            = sva_res$dds,
      comparisons        = comparison_results_list
    )
  }


# fit <- lmFit(logCPM, design)
# fit <- treat(fit, lfc=log2(1.2), trend=TRUE)
# topTreat(fit, coef=ncol(design))



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
remove_outliers.DGEList <-
  function(
    object,
    design,
    pc1_zscore_cutoff,
    pc2_zscore_cutoff = NULL
  ){

    v <- voom(object, design)
    pca_res =
      prcomp_irlba(v$E)[['rotation']] %>%
      as_tibble() %>%
      mutate(
        sample = colnames(v$E),
        zscore = abs((PC1 - mean(PC1))/sd(PC1))
      ) %>%
      inner_join(as_tibble(v$targets, rownames = "sample")) %>%
      rename(sample_name = sample)

    pc1_outliers <-
      pca_res %>%
      filter(PC1 >= pc1_zscore_cutoff) %>%
      pull(sample_name)

    if (!is.null(pc2_zscore_cutoff)){
      pc2_outliers <-
        pca_res %>%
        filter(PC2 >= pc2_zscore_cutoff) %>%
        pull(sample_name)
    } else {
      pc2_outliers <- NULL
    }

    outliers <- unique(c(pc1_outliers, pc2_outliers))

    if (length(outliers > 0)){
      object <- object[,colnames(object) %nin% outliers]
    }

    return(
      list(
        count_object = object,
        pca          = pca_res,
        removed      = outliers
      )
    )
  }


#' @rdname remove_outliers
#' @method remove_outliers DESeqDataSet
#' @importFrom irlba prcomp_irlba
#' @importFrom DESeq2 vst estimateSizeFactors
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate inner_join filter pull
#' @return DESeqDataSet
remove_outliers.DESeqDataSet <-
  function(
    object,
    pc1_zscore_cutoff,
    pc2_zscore_cutoff = NULL
  ){

    object <-
      estimateSizeFactors(
        object,
        locfun = shorth,
        type = "poscounts")

    vsd <- assay(vst(object))
    pca_res = prcomp_irlba(vsd)[['rotation']] %>%
      as_tibble() %>%
      mutate(
        sample = colnames(vsd),
        zscore = abs((PC1 - mean(PC1))/sd(PC1))) %>%
      inner_join(as_tibble(colData(object), rownames = "sample"))

    pc1_outliers <- pca_res %>%
      filter(pc1_zscore >= pc1_zscore_cutoff) %>%
      pull(sample_name)

    if (!is.null(pc2_zscore_cutoff)){
      pc2_outliers <- pca_res %>%
        filter(pc2_zscore >= pc2_zscore_cutoff) %>%
        pull(sample_name)
    } else {
      pc2_outliers <- NULL
    }

    outliers <- unique(c(pc1_outliers, pc2_outliers))

    if (length(outliers > 0)){
      dds <- dds[,colnames(dds) %nin% outliers]
    }

    return(
      list(
        count_object = object,
        pca          = pca_res,
        removed      = outliers
      )
    )
  }



###--- kmeans version ---###
sample_clustering <- function(
  exprs_mat,
  from = 2,
  to = 20,
  by = 1){

  kmeans_clusters <-
    future_map(seq(from = from,
                   to = to,
                   by = by),
               function(i){
                 clara(
                   x = exprs_mat,
                   k = i,
                   metric = "jaccard",
                   stand = TRUE,
                   samples = 50,
                   pamLike = TRUE
                 )
               }
    ) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by))

  sils <- map(kmeans_clusters, function(j){
    silhouette(j)
  })

  avg_sils <- map_dbl(sils, function(k){
    mean(k[,3])
  }) %>%
    set_names(seq(from = from,
                  to = to,
                  by = by)) %>%
    enframe() %>%
    mutate(name = as.integer(name))

  optimal_k =
    avg_sils %>%
    filter(name > 2) %>%
    top_n(
      n = 1,
      wt = value
    )

  clusters <- kmeans_clusters[[as.character(optimal_k[["name"]])]][["clustering"]] %>%
    enframe(name = "sample_name",
            value = "cluster") %>%
    mutate(cluster = as_factor(cluster))

  return(list(avg_sils = avg_sils,
              optimal_k = optimal_k,
              clusters = clusters))
}

# Use the output from cluster_silhouette()
plot_resolution_silhouette_coeff <- function(cluster_silhouette){
  cluster_silhouette %>%
    ggplot(aes(x = res,
               y = coeff)) +
    geom_line() +
    theme_cowplot()
}

ident_clusters <- function(expr_mat,
                           optimal_k_method = "Tibs2001SEmax",
                           nstart = 25,
                           K.max = 50,
                           B = 100,
                           d.power = 2){

  message("Performing random forest")
  module_rf <-
    randomForest(
      x = expr_mat,
      y = NULL,
      prox = T)

  message("Calculating sample distances from random forest proximity")
  rf_distance_mat <-
    parallelDist::parallelDist(1 - module_rf$proximity) %>%
    as.matrix()

  message("Calculating goodness of k-means clustering")
  kmeans_gap_stat <-
    clusGap(
      x = rf_distance_mat,
      FUNcluster = kmeans,
      nstart = nstart,
      K.max = K.max,
      B = B,
      d.power = d.power,
      parallel = TRUE,
      future_plan = "multisession"
    )

  message("Optimizing k")
  new_optimal_k <-
    with(
      data = kmeans_gap_stat,
      expr = maxSE(Tab[,"gap"],
                   Tab[,"SE.sim"],
                   method=optimal_k_method
      )
    )

  message("Calculating clusters")
  k_clusters <-
    kmeans(
      x = rf_distance_mat,
      centers = new_optimal_k,
      nstart = 25
    )

  sample_clusters <-
    enframe(x = k_clusters[["cluster"]],
            name = "sample_name",
            value = "cluster")

  list(
    kmeans_res = k_clusters,
    rf_distance = rf_distance_mat,
    clusters = sample_clusters,
    gap_stat = kmeans_gap_stat
  )
}

plot_dispersion_estimate <- function(object, CV = FALSE){
  px <- mcols(object)$baseMean
  sel <- (px > 0)
  px <- px[sel]
  f <- ifelse(CV, sqrt, I)
  py <- f(mcols(object)$dispGeneEst[sel])
  ymin <- 10^floor(log10(min(py[py > 0], na.rm = TRUE)) -
                     0.1)

  outlier_shape <- ifelse(mcols(object)$dispOutlier[sel],
                          1, 16)
  outlier_size <- ifelse(mcols(object)$dispOutlier[sel],
                         2 * 0.45, 0.45)
  outlier_halo <- ifelse(mcols(object)$dispOutlier[sel],
                         "final", "gene-est")

  disp_data <- tibble(px = px,
                      py = pmax(py, ymin),
                      outlier_shape = as_factor(outlier_shape),
                      outlier_size = as_factor(outlier_size),
                      outlier_halo = as_factor(outlier_halo),
                      dispersions = f(dispersions(object)[sel]),
                      dispersions_fit = f(mcols(object)$dispFit[sel]))

  disp_plot <- disp_data %>%
    ggplot(aes(x = px,
               y = py)) +
    geom_point() +
    geom_point(aes(x = px,
                   y = dispersions,
                   size = outlier_size,
                   shape = outlier_shape,
                   color = outlier_halo)) +
    scale_x_log10() +
    scale_y_log10() +
    scale_shape_manual(values = c(1, 16)) +
    scale_size_manual(values = c(1,2)) +
    scale_color_manual(values = c(
      "dodgerblue",
      "red",
      "black"), ) +
    geom_line(mapping = aes(x = px,
                            y = dispersions_fit,
                            color = "fitted"),
              size = 1) +
    labs(x = "mean of normalized counts",
         y = "dispersion",
         color = "") +
    guides(size = "none",
           shape = "none") +
    theme_cowplot() +
    theme(legend.justification=c(1,0), legend.position=c(1,0))

  return(disp_plot)
}

fix_antibody_values <- function(i) {
  recode(.x = i,
         Negative = "negative",
         POSITIVE = "positive",
         `no_val` = "no_val",
         Indeterminate = "indeterminate")
}


#' Title
#'
#' @param object DESeqResults object
#'
#' @return
#' @export
#'
#' @examples
alt_summary <- function(object){
  notallzero <- sum(object$baseMean > 0)
  up <- sum(object[["padj"]] < 0.05 & object$log2FoldChange >
              metadata(object)$lfcThreshold, na.rm = TRUE)
  down <- sum(object[["padj"]] < 0.05 & object$log2FoldChange <
                metadata(object)$lfcThreshold, na.rm = TRUE)
  outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
  if (is.null(metadata(object)$filterThreshold)) {
    ft <- 0
  } else {
    ft <- round(metadata(object)$filterThreshold)
  }

  filt <- sum(!is.na(object$pvalue) & is.na(object$padj))

  total <- nrow(object)


  tibble(
    up = up,
    down = down,
    outlier = outlier,
    ft = ft,
    lowcounts = filt,
    total = total
  )
}

calc_sva <- function(object, ...){
  UseMethod("calc_sva")
}

calc_sva.DESeqDataSet <- function(object, model_design = NULL, n.sva = NULL){
  model_design_factors <-
    model_design %||% as.character(design(object))[[2]] %>%
    str_remove(pattern = "~")

  n.sva <- n.sva %||% 2

  dat <- counts(object, normalized = TRUE) %>%
    purrr::keep(.p = ~ rowMeans(.x) > 1)

  model_design <-
    as.formula(
      paste("~",
            paste(
              unlist(model_design_factors),
              collapse = " + "),
            collapse = " ")
    )

  mod  <- model.matrix(design(object), colData(object))

  mod0 <- model.matrix(~ 1, colData(object))

  svseq <- svaseq(dat, mod, mod0, n.sv = n.sva)

  colnames(svseq$sv) <- paste0("SV", seq(ncol(svseq$sv)))

  for (i in seq(ncol(svseq$sv))){
    object[[paste0("SV",i)]] <- svseq$sv[,i]
  }

  design(object) <-
    as.formula(
      paste("~",
            paste(model_design_factors,
                  paste(colnames(svseq$sv),
                        collapse = " + "),
                  sep = " + "),
            collapse = " "))

  object <- DESeq(object, parallel = TRUE)

  ret_vals = list(
    dds = object,
    sva = svseq
  )
}

calc_sva.DGEList <- function(object, model_design = NULL, batch_var = NULL, n.sva = 2){

  batch_var <- batch_var %||% 1

  svseq <-
    svaseq(
      dat = object[["counts"]],
      mod = model.matrix(as.formula(model_design), object[["samples"]]),
      mod0 = model.matrix(as.formula(paste("~", 1)), object[["samples"]]),
      n.sv = n.sva
    )

  svseq[["sv"]] <-
    set_names(
      x = as_tibble(svseq[["sv"]], .name_repair = "unique"),
      nm = paste0("SV", seq(ncol(svseq[["sv"]])))
    ) %>%
    mutate(sample_name = rownames(object[["samples"]]))

  object[["samples"]] <-
    left_join(
      as_tibble(object[["samples"]], rownames = "sample_name"),
      svseq[["sv"]]
    ) %>%
    column_to_rownames(var = "sample_name")

  svseq$sv <- column_to_rownames(svseq[["sv"]],"sample_name")

  design_formula <-
    as.formula(
      paste("~",
            paste(as.character(model_design)[[2]],
                  paste(colnames(svseq[["sv"]]),
                        collapse = " + "),
                  sep = " + "),
            collapse = " "))

  ret_vals = list(
    dge = object,
    sva = svseq,
    design = design_formula
  )
}

plot_sva <- function(sva_graph_data){
  sva_graph_data %>%
    as_tibble(rownames = "sample_name") %>%
    select(
      sample_name,
      starts_with("SV")
    ) %>%
    pivot_longer(
      -sample_name,
      names_to = "covar"
    ) %>%
    ggplot(
      aes(
        x = sample_name,
        y = value
      )
    ) +
    geom_point() +
    geom_hline(
      yintercept = 0,
      color = "red"
    ) +
    facet_grid(
      rows = vars(covar)
    ) +
    theme_cowplot() +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          size = 9,
          hjust = 1,
          vjust = 1
        )
    )
}

drake_recode <- function(target_list, thing_to_unquote_splice){
  dplyr::recode(.x = target_list, !!! {{thing_to_unquote_splice}})
}

# An improved version of rstatix::add_y_position
# doesn't take 30 minutes to run either
grouped_add_xy_positions <- function(stats_tbl,
                                     data_tbl,
                                     group_var,
                                     compare_value,
                                     cutoff = 0.05,
                                     step_increase = 0.1){

  unique_groups <- stats_tbl %>% pull({{group_var}}) %>% unique()

  data_min_max <-
    data_tbl %>%
    select({{group_var}}, {{compare_value}}) %>%
    group_by({{group_var}}) %>%
    summarise(max = max({{compare_value}}),
              min = min({{compare_value}}),
              span = max-min,
              step = span * step_increase)

  tbl_with_positions <- map_dfr(unique_groups, function(x){
    stats_subset <- stats_tbl %>% filter({{group_var}} == x) %>% add_x_position()

    stats_subset <- if ("p.adj" %in% names(stats_subset)){
      stats_subset %>% filter(p.adj <= cutoff)
    } else {
      stats_subset %>% filter(p <= cutoff)
    }

    min_max_subset <- data_min_max %>% filter({{group_var}} == x)
    if (nrow(stats_subset) > 1){
      positions <-
        seq(
          from = min_max_subset[['max']],
          by = min_max_subset[['step']],
          to = min_max_subset[['max']] + nrow(stats_subset)*min_max_subset[['step']])
      stats_subset[['y.position']] <- positions[2:length(positions)]
    }
    stats_subset
  })
  return(tbl_with_positions)
}

convert_nuID_to_probeID <- function(object, rename_list){
  featureNames(object) <-
    featureNames(object) %>%
    recode(!!!rename_list)

  return(object)
}


#' @rdname add_sample_metadata
#' @export
add_sample_metadata <- function(object, ...)
  UseMethod(generic = "add_sample_metadata")

#' @rdname add_sample_metadata
#' @alias add_sample_metadata,DGEList-method
#'
#' @param object a \code{DGEList} object
#' @param md a \code{data.frame} with one or more columns containing data to add.
#' MUST have the same rownames as the samples in the DGEList
#'
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble column_to_rownames
#'
#' @exportMethod "add_module_score"
add_sample_metadata.DGEList <-
  function(object, md){
    object$samples <-
      left_join(
        x = as_tibble(object$samples, rownames = "rownames"),
        y = as_tibble(md, rownames = "rownames")
      ) %>%
      column_to_rownames(var = "rownames")
    object
  }

#' @rdname add_sample_metadata
#' @alias add_sample_metadata,DESeqDataSet-method
#'
#' @param object a \code{DESeqDataSet} object
#' @param md a \code{data.frame} with one or more columns containing data to add.
#' MUST have the same rownames as the samples in the DESeqDataSet
#'
#' @importFrom dplyr left_join
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors DataFrame
#'
#' @exportMethod "add_module_score"
add_sample_metadata.DESeqDataSet <-
  function(object, scores) {
    colData(object) <-
      left_join(
        x = as_tibble(colData(object), rownames = "rownames"),
        y = as_tibble(scores, rownames = "rownames")
      ) %>%
      DataFrame()

    object
  }


#' @rdname retrieve_metadata
#' @export
retrieve_metadata <- function(object, ...){
  UseMethod("retrieve_metadata")
}

#' @rdname retrieve_metadata
#' @alias retrieve_metadata,DGEList-method
#'
#' @param object a \code{DGEList} object
#' @param columns a character list or numeric indices of the columns to retrieve.
#' If NULL, all columns are retrieved.  Default: NULL
#'
#' @importFrom rlang %||%
#' @importFrom tibble as_tibble
#'
#' @exportMethod "retrieve_metadata"
retrieve_metadata.DGEList <- function(object, columns = NULL){
  metadata <- as_tibble(object$samples, rownames = "sample_names")
  columns <- columns %||% colnames(metadata)
  metadata[,columns]
}

#' @rdname retrieve_metadata
#' @alias retrieve_metadata,DESeqDataSet-method
#'
#' @param object a \code{DESeqDataSet} object
#' @param columns a character list or numeric indices of the columns to retrieve.
#' If NULL, all columns are retrieved.  Default: NULL
#'
#' @importFrom rlang %||%
#' @importFrom tibble as_tibble
#'
#' @exportMethod "retrieve_metadata"
retrieve_metadata.DESeqDataSet <- function(object, columns = NULL){
  metadata <- as_tibble(colData(object), rownames = "sample_names")
  columns <- columns %||% colnames(metadata)
  metadata[,columns]
}


# Adapted from Seurat::AddModuleScore, which in turn took it from Tirosh (2006)
tirosh_score_modules <- function(expr_obj, module_list, breaks = 25, num_ctrls = 100) {

  features       <- module_list
  name           <- "module"
  cluster_length <- length(x = features)

  data_avg       <- Matrix::rowMeans(x = expr_obj)
  data_avg       <- data_avg[order(data_avg)]
  data_cut       <-
    cut_number(
      x      = data_avg + rnorm(n = length(data_avg)) / 1e30,
      n      = num_ctrls,
      labels = FALSE,
      right  = FALSE
    )

  names(x = data_cut) <- names(x = data_avg)
  ctrl_use <-
    vector(
      mode   = "list",
      length = cluster_length
      )

  # for each module
  for (i in seq(cluster_length)) {
    # use only the module genes that are present in our dataset
    features_use <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]

    # for each module gene
    for (j in seq_along(features_use)) {
      ctrl_use[[i]] <-
        c(
          ctrl_use[[i]],
          names(
            x = sample(
              x       = data_cut[which(x = data_cut == data_cut[features_use[j]])],
              size    = num_ctrls,
              replace = FALSE
            )
          )
        )
    }
  }

  ctrl_use    <- lapply(X = ctrl_use, FUN = unique)
  ctrl_scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl_use),
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:length(ctrl_use)) {
    features_use     <- ctrl_use[[i]]
    ctrl_scores[i, ] <- Matrix::colMeans(x = expr_obj[features_use, ])
  }

  features_scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster_length,
    ncol = ncol(x = expr_obj)
  )

  for (i in 1:cluster_length) {
    features_use         <- features[[i]][which(features[[i]] %in% rownames(expr_obj))]
    data_use             <- expr_obj[features_use, , drop = FALSE]
    features_scores[i, ] <- Matrix::colMeans(x = data_use)
  }

  features_scores_use <- features_scores - ctrl_scores
  rownames(x = features_scores_use) <- names(module_list)
  features_scores_use <- as.data.frame(x = t(x = features_scores_use))
  rownames(x = features_scores_use) <- colnames(x = expr_obj)

  features_scores_use
}

score_modules <- function(res, modules){
  res %>%
    as_tibble(rownames = "gene") %>%
    inner_join(gene_module_tbl) %>%
    filter(padj < 0.05) %>%
    group_by(module) %>%
    summarise(
      percent_pos = sum(log2FoldChange > 0) / n(),
      percent_neg = sum(log2FoldChange < 0) / n(),
      percent_diff = percent_pos - percent_neg) %>%
    dplyr::select(module, percent_diff) %>%
    rename({{compare_class}} := percent_diff)
}


#' @title find_sequencing_samples
#'
#' @param seq_file_directory
#' @param source_names
#' @param alignment_source
#'
#' @return
#' @export
#'
#' @importFrom purrr set_names discard pluck map_chr
#' @importFrom stringr str_remove str_split
#' @importFrom janitor make_clean_names
#'
#' @examples
find_sequencing_samples <- function(
  seq_file_directory,
  source_names = NULL,
  alignment_source = "salmon"){

  tx_files <-
    switch(
      EXPR = alignment_source,
      salmon = find_sequencing_samples.salmon(seq_file_directory = seq_file_directory, source_names = source_names)
    )

  return(tx_files)
}

find_sequencing_samples.salmon <-
  function(
    seq_file_directory,
    source_names = NULL){
    tx_files <-
      dir(
        path = seq_file_directory,
        pattern = "quant.sf.gz",
        recursive = TRUE,
        full.name = TRUE
      ) %>%
      grep(
        pattern = "Undetermined|NONE",
        invert = TRUE,
        value = TRUE
      ) %>%
      set_names(
        x = .,
        nm = . %>%
          str_split(pattern = "/") %>%
          map_chr(~pluck(.x, length(.x)-1)) %>%
          str_remove(pattern = '(_[L|S][[:digit:]]+)+') %>%
          janitor::make_clean_names(case = "all_caps")
      )

    if (!is.null(source_names)){
      tx_files <-
        purrr::discard(
          .x = tx_files,
          .p =
            is.na(
              match(
                names(tx_files),
                source_names
              )
            )
        )
    }
    tx_files
  }
