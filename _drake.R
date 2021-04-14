options(future.globals.maxSize = +Inf)
source("R/packages.R")        # Load your packages, e.g. library(drake).
source("R/functions.R")       # Define your custom code as a bunch of functions.
source("R/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules

c5 <- read.gmt("references/c5.all.v6.2.symbols.gmt")
Sys.setenv('RSTUDIO_PANDOC' = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

load(file="references/gencode_v32_virus_tx2gene_v1.2.RData") #generated from rtracklayer::readGFF()

### Setup project variables
projects_to_include             = NULL # "BChong2019.1"
projects_to_exclude             = "BChong2019.1"
disease_classes_to_include      = c("Control", "SLE")
disease_classes_to_exclude      = NULL
study_design                    = ~ run_id + disease_class
comparison_grouping_variable    = "disease_class"
control_group                   = "Control"
experimental_group              = "SLE"
batch_variable                  = "run_id"
num_sva                         = 3
analysis_library                = "limma"

initial_concentration_threshold = 1.5
pc1_zscore_threshold            = 2
pc2_zscore_threshold            = 2

### Setup file locations
seq_file_directory              = "/home/rstudio/workspace/datasets/rnaseq/novaseq/"
metadata_file                   = "metadata/NovaSeq_Sample_List.xlsx"
metadata_sheet                  = "Sheet1"
clinical_file                   = ""

BPPARAM = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(strategy = future::multisession)

analysis_plan = code_to_plan("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(
  plan                          = analysis_plan,
  verbose                       = 2,
  parallelism                   = "future",
  log_progress                  = TRUE,
  jobs                          = parallel::detectCores()-2,
  lock_envir                    = FALSE
  )
