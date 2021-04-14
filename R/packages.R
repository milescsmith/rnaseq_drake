reticulate::use_condaenv('reticulate', required = TRUE, conda = "~/conda/bin/conda")
required_packages <-
  c(
    "caret",
    "clusterProfiler",
    "corrplot",
    "cowplot",
    "data.table",
    "DESeq2",
    "drake",
    "edgeR",
    "factoextra",
    "flextable",
    "formattable",
    "furrr",
    "ggforce",
    "ggplotify",
    "ggpubr",
    "ggrepel",
    "ggtext",
    "gtools",
    "HGNChelper",
    "irlba",
    "janitor",
    "limma",
    "kableExtra",
    "knitr",
    "magrittr",
    "matrixStats",
    "moduleScoreR",
    "paletteer",
    "parallelCluster",
    "pheatmap",
    "plotly",
    "RColorBrewer",
    "readxl",
    "rlang",
    "rstatix",
    "scales",
    "stats",
    "sva",
    "tidymodels",
    "randomForest",
    "tidyverse",
    "tximport",
    "uwot",
    "viridis",
    "WGCNA"
  )

lapply(required_packages, function(x){
  if (x %in% installed.packages()){
    require(x, character.only = TRUE)
  } else {
    BiocManager::install(
      pkgs = x,
      Ncpus = parallel::detectCores(),
      dependencies = c("Depends","Imports","LinkingTo"),
      clean = TRUE,
      upgrade = FALSE,
      ask = FALSE,
      checkBuilt = TRUE
      )
    require(x, character.only = TRUE)
  }
})

