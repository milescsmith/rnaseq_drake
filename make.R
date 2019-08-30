source("R/packages.R")  # Load your packages, e.g. library(drake).
source("R/functions.R") # Define your custom code as a bunch of functions.


### Setup bucket access
flyio_set_datasource("gcs")
flyio_auth("/opt/google_project_scrna_196615_storage_key.json")
flyio_set_bucket("memory_alpha", data_source="gcs")

import_rda(file="references/gencode.v31.annotation.gtf.RData",
           bucket = "memory_alpha") #generated from rtracklayer::readGFF()

import_rda(file="references/banchereau_modules.RData",
           bucket = "memory_alpha",
           data_source = "gcs")

### Setup project variables
projects_to_include = c("Xencor")
projects_to_exclude = c("none")
study_design = ~ rin + initial_concentration_ng_ul + run_id + disease_class
comparison_grouping_variable = "disease_class"
control_group = "Control"
experimental_group = "SLE"

### Setup file locations
seq_file_directory = "/home/milo/datasets/S4"
metadata_file = "datasets/rnaseq/S4/NovaSeq_Sample_List.xlsx"

source("R/plan.R")      # Create your drake plan.
config <- drake_config(plan)
vis_drake_graph(config)

BPPARAM = BiocParallel::MulticoreParam(workers=parallel::detectCores())
BiocParallel::register(BPPARAM)

make(plan = plan,
     verbose = 2,
     parallelism = "future",
     jobs = 8,
     lock_envir = FALSE)
