# This file serves the r_*() functions (e.g. r_make()) documented at
# https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity # nolint
# and
# https://ropensci.github.io/drake/reference/r_make.html

# Load your packages and supporting functions into your session.
# If you use supporting scripts like the ones below,
# you will need to supply them yourself. Examples:
# https://github.com/wlandau/drake-examples/tree/master/main/R
source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.

### Setup bucket access
flyio_set_datasource("gcs")
flyio_auth("/opt/google_project_scrna_196615_storage_key.json")
flyio_set_bucket("memory_alpha", data_source="gcs")

import_rda(file="references/gencode.v31_viruses_tx2gene.RData",
           bucket = "memory_alpha") #generated from rtracklayer::readGFF()

import_rda(file="references/banchereau_modules.RData",
           bucket = "memory_alpha",
           data_source = "gcs")

import_rda(file="references/kegerreis_ldg_modules.RData",
           bucket = "memory_alpha",
           data_source = "gcs")

### Setup project variables
projects_to_include = c("ABC")
projects_to_exclude = c("Xencor")
disease_classes_to_include = c("Control", "SLE")
disease_classes_to_exclude = NULL
study_design = ~ initial_concentration_ng_ul + run_id + disease_class
comparison_grouping_variable = "disease_class"
control_group = "Control"
experimental_group = "SLE"
initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2

### Setup file locations
seq_file_directory = "/home/milo/charon/datasets/novaseq"
metadata_file = "datasets/rnaseq/novaseq/NovaSeq_Sample_List.xlsx"

source("R/plan.R")
# config <- drake_config(plan)
# vis_drake_graph(config)

BPPARAM = SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)

# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan,
             verbose = 2,
             parallelism = "future",
             jobs = parallel::detectCores(),
             lock_envir = TRUE)
