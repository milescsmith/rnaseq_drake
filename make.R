source("/home/milo/datasets/bulk_preclinical/rnaseq_drake_salmon/R/packages.R")
source("/home/milo/datasets/bulk_preclinical/rnaseq_drake_salmon/R/functions.R")


### Setup bucket access
flyio_set_datasource("gcs")
flyio_auth("/opt/google_project_scrna_196615_storage_key.json")
flyio_set_bucket("memory_alpha", data_source="gcs")

import_rda(file="references/gencode.v31_viruses_tx2gene.RData",
           bucket = "memory_alpha") #generated from rtracklayer::readGFF()

import_rda(file="references/banchereau_modules.RData",
           bucket = "memory_alpha",
           data_source = "gcs")

### Setup project variables
study_design = ~ ethnicity + celltype + disease_class
comparison_grouping_variable = "disease_class"

### Setup file locations
seq_file_directory = "/home/milo/datasets/bulk_preclinical/salmon"
sgl <- read_csv("~/datasets/bulk_preclinical/sams_rnaseq_genelist.csv", col_names = c("gene","pathway"))
sfgl <- read_xlsx("~/datasets/bulk_preclinical/sams_full_rnaseq_genelist.csv.xlsx")

source("/home/milo/datasets/bulk_preclinical/rnaseq_drake_salmon/R/plan.R")      # Create your drake plan.
# config <- drake_config(plan)
# vis_drake_graph(config)

BPPARAM = BiocParallel::MulticoreParam(workers=parallel::detectCores())
BiocParallel::register(BPPARAM)

make(plan = plan,
     verbose = 2,
     parallelism = "future",
     jobs = 8,
     lock_envir = FALSE)
