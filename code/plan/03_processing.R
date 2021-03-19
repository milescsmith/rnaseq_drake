vsd = vst(dds_processed)

vsd_exprs = assay(vsd)

top_vars =
  rowSds(vsd_exprs) %>%
  set_names(rownames(vsd_exprs)) %>%
  enframe() %>%
  top_n(20000, value) %>%
  pull(name)

vsd_top = vsd_exprs[top_vars,] %>%
  t()
