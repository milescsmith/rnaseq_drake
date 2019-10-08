#!/usr/bin/env Rscript

BPPARAM = BiocParallel::SnowParam(workers=96, type="SOCK")
BiocParallel::register(BPPARAM)
library(drake)
r_make()