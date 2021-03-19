sample_dists =
  vsd_exprs %>%
  t() %>%
  parallelDist::parallelDist()

#fig.width=12, fig.height=9
sampleDistMatrix = as.matrix(sample_dists)

sample_dendrogram =
  sample_dists %>%
  hclust() %>%
  as.dendrogram()

# okay, so we *were* running this on the annotated modules
# but that seems short sighted.  Ignores genes that might fall into
# an unannotated module or just nothing.
# changing to running on the top 20000 variable genes instead.
sample_cluster_info =
  ident_clusters(
    vsd_exprs,
    K.max = 20
  )

# clusters =
#   leiden_cluster(
#     exprs = vsd_top,
#     res = 1.5,
#     nneighbors = 10,
#     column_name = "cluster"
#   )

clusters =
  sample_cluster_info$clusters %>% mutate(cluster = as_factor(cluster))

study_md =
  colData(dds_with_scores) %>%
  as_tibble(rownames = "sample_name") %>%
  left_join(clusters)

annotation_info =
  select(.data = study_md,
         disease_class,
         sex,
         cluster,
         sample_name) %>%
  column_to_rownames(var = "sample_name")
