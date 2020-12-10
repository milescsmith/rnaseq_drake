options(future.globals.maxSize = +Inf)
source("R/packages.R")        # Load your packages, e.g. library(drake).
source("R/functions.R")       # Define your custom code as a bunch of functions.
source("R/extant_modules.R")  # Banchereau, Kegerreis, and Metagene modules

c5 <- read.gmt("references/c5.all.v6.2.symbols.gmt")
Sys.setenv('RSTUDIO_PANDOC' = '/usr/lib/rstudio-server/bin/pandoc')
WGCNA::allowWGCNAThreads()

load(file="references/gencode_v32_virus_tx2gene_v1.2.RData") #generated from rtracklayer::readGFF()


### Setup project variables
study_design = ~ cell_type + disease_class
comparison_grouping_variable = "disease_class"
control_group = "control"
experimental_group = "SLE"
batch_variable = "plate_id"
num_sva = 3

initial_concentration_threshold = 1.5
pc1_zscore_threshold = 2
pc2_zscore_threshold = 2

### Setup file locations
seq_file_directory = "/home/rstudio/workspace/datasets/rnaseq/narch_advanta/data/results/salmon"
metadata_file = "metadata/NARCH_deidentification.xlsx"
metadata_sheet = "Sheet1"

BPPARAM = BiocParallel::SnowParam(workers=parallel::detectCores(), type = "SOCK")
BiocParallel::register(BPPARAM)
future::plan(future::multisession)

source("R/plan.R")
# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake_config(plan = analysis_plan,
             verbose = 2,
             parallelism = "future",
             log_progress = TRUE,
             jobs = parallel::detectCores()-2,
             lock_envir = FALSE)
