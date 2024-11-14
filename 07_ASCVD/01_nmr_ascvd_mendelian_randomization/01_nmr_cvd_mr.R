#!/usr/bin/env Rscript

## script to run MR of NMR traits on major CVD outcomes
## This script is handled by SLURM via `03_nmr_cvd_mr_submit.sh`

## Carl Beuchel 2024-08-06
rm(list = ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly = T)

## little options
options(stringsAsFactors = F)

# SETUP ----

# Directory setup for saving results and graphics
setwd("</YOUR/PROJECT/DIRECTORY>")

## --> packages required <-- ##
for (i in c(
  "data.table",
  "magrittr",
  "here",
  "data.table",
  "glue",
  "doMC",
  "TwoSampleMR",
  "cowplot",
  "fs"
)
) {
  suppressPackageStartupMessages(
    library(
      i, 
      character.only = TRUE
    ))
}

plot_save_path <- here("graphics/02_nmr_cvd_mr/")
result_save_path <- here("output/02_nmr_cvd_mr/")
dir.create(plot_save_path, showWarnings = FALSE)
dir.create(result_save_path, showWarnings = FALSE)

# INPUT ----

# metabolite_index <- 180L # for DEBUG
metabolite_index <- as.integer(args[1])

if ((length(metabolite_index) != 1) & (!is.integer(metabolite_index))) {
  stop("Something went wrong with the index!")
}

# NMR GWAS Summary Statistics
exposure_files <- fs::dir_ls("03_nmr_cvd_mr/input/exposure_data/")

# ASCVD Summary Statistics
outcome_files <- fs::dir_ls("03_nmr_cvd_mr/input/outcome_data/")

# Different Variant Sets trimmed for pleiotropic instruments
variant_set_files <- fs::dir_ls("03_nmr_cvd_mr/input/variant_sets/")

# Annotation File containing meta data for each ASCVD outcome
outcome_annotation <- fread("/PATH/TO/ANNOTATION/FILE.csv", na.strings = "")

# MAIN ----

## EXPOSURE DATA ----

# Only One Metabolite
exposure_file <- exposure_files[metabolite_index]
exposure_id <- gsub("_allchr.txt.gz", "", exposure_file, fixed = TRUE)
exposure_id <- gsub(
  "03_nmr_cvd_mr/input/exposure_data/gwas_", 
  "",
  exposure_id, 
  fixed = TRUE
)
exposure_data <- read_exposure_data(
  filename = exposure_file,
  clump = FALSE,
  snp_col = "ID",
  chr_col = "CHROM",
  pos_col = "GENPOS",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "LOG10P",
  eaf_col = "A1FREQ",
  other_allele_col = "ALLELE0",
  effect_allele_col = "ALLELE1",
  log_pval = TRUE,
  samplesize_col = "N"
)
setDT(exposure_data)
exposure_data[, chr.exposure := as.character(chr.exposure)]
exposure_data[chr.exposure == "23", chr.exposure := "X"]
exposure_data[, ID := SNP]
exposure_data[
  ,
  SNP := paste(
    chr.exposure, 
    pos.exposure, 
    pmin(other_allele.exposure, effect_allele.exposure),
    pmax(other_allele.exposure, effect_allele.exposure), 
    sep = "_"
  )
]

## MAIN LOOP ----

# Test across outcomes and variant sets
test_grid <- expand.grid(
  outcome = outcome_files,
  variant_set = variant_set_files, 
  stringsAsFactors = FALSE
)
setDT(test_grid)
test_grid[, index := .I]

cl <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)
all_outcomes_mr_results <- foreach(
  i = test_grid$index
) %dopar% {
  
  # Set exposure and outcome ID
  outcome_file <- test_grid[i, outcome]
  variant_set_file <- test_grid[i, variant_set]
  
  # Select variant set
  variant_set_data <- fread(variant_set_file)
  mr_variants <- variant_set_data[exposure == exposure_id, marker_name]
  variant_set_id <- variant_set_data$variant_set[1]
  
  # outcome data
  outcome_id <- gsub("_grch37.tsv.bgz", "", outcome_file, fixed = TRUE)
  outcome_id <- gsub(
    "/?03_nmr_cvd_mr/input/outcome_data/+", 
    "", 
    outcome_id
  )
  outcome_data <- read_outcome_data(
    filename = outcome_file,
    sep = "\t",
    snp_col = 'SNP',
    chr_col = 'CHR',
    pos_col = 'BP',
    other_allele_col = 'A1',
    effect_allele_col = 'A2',
    beta_col = 'BETA',
    se_col = 'SE',
    pval_col = 'P',
    eaf_col = 'FRQ',
    log_pval = FALSE
  )
  setDT(outcome_data)
  outcome_data[, ID := SNP]
  outcome_data[
    ,
    SNP := paste(
      chr.outcome, 
      pos.outcome, 
      pmin(other_allele.outcome, effect_allele.outcome),
      pmax(other_allele.outcome, effect_allele.outcome), 
      sep = "_"
    )
  ]
  
  # Filter data to SNP in variant set
  filtered_exposure_data <- exposure_data[SNP %in% mr_variants] %>% copy()
  outcome_data <- outcome_data[SNP %in% mr_variants]
  mr_data <- harmonise_data(
    exposure_dat = filtered_exposure_data, 
    outcome_dat = outcome_data, 
    action = 2
  )
  
  # Check if one of the datasets is empty
  if (nrow(mr_data) <= 10) {
    message("[LOG] Not enough overlapping SNP to continue.")
    return(data.table())
  }
  
  # Annotate sample size for outcome data
  mr_data$samplesize.outcome <- outcome_annotation[phenotype_id == outcome_id, sample_size]
  mr_data$ncase.outcome <- outcome_annotation[phenotype_id == outcome_id, n_cases]
  mr_data$ncontrol.outcome <- outcome_annotation[phenotype_id == outcome_id, n_controls]
  
  ## MR ----
  
  mr_data <- steiger_filtering(mr_data)
  setDT(mr_data)
  
  # H0 is association with outcome rather than exposure, so keep sig. results
  mr_data <- mr_data[steiger_pval < 0.05, ]
  
  if (nrow(mr_data) <= 10) {
    message("[LOG] Not enough overlapping SNP to continue.")
    return(data.table())
  }
  
  setDF(mr_data)
  mr_results <- mr(
    mr_data, 
    method_list = c("mr_egger_regression", "mr_ivw")
  )
  
  # Test for pleiotropy (MR-Egger Intercept)
  mr_results_pleiotropy <- mr_pleiotropy_test(mr_data)
  
  # Calcualate the heterogeneity
  mr_results_heterogeneity <- mr_heterogeneity(mr_data)
  
  ## PLOT ----
  
  # Scatter Plot
  p1 <- mr_scatter_plot(mr_results, mr_data)
  
  # Forest Plot
  res_single <- mr_singlesnp(mr_data)
  p2 <- mr_forest_plot(res_single)
  
  # Funnel Plot
  p3 <- mr_funnel_plot(res_single)
  
  res_plot <- plot_grid(
    plotlist = list(p1[[1]], p3[[1]], p2[[1]]),
    ncol = 3,
    nrow = 1
  )
  plot_filename <- glue(
    "{plot_save_path}/mr_plots_{exposure_id}_{outcome_id}_{variant_set_id}.png"
  )
  ggsave2(
    filename = plot_filename,
    plot = res_plot, 
    width = 10, 
    height = 4,
  )
  
  ## RETURN ----
  
  # Add useful information to this
  setDT(mr_results)
  setDT(mr_results_heterogeneity)
  setDT(mr_results_pleiotropy)
  mr_results_heterogeneity$outcome <- NULL
  mr_results_heterogeneity$exposure <- NULL
  mr_results$outcome <- NULL
  mr_results$exposure <- NULL
  
  # Intercept is estimate for pleiotropy
  mr_results[
    method == "MR Egger", 
    c("b_intercept", "se_intercept", "p_intercept") := 
      mr_results_pleiotropy[,.(egger_intercept, se, pval)]
    ]
  
  # Add heterogeneity test results
  setkey(mr_results_heterogeneity, id.exposure, id.outcome, method)
  setkey(mr_results, id.exposure, id.outcome, method)
  mr_results <- mr_results[mr_results_heterogeneity]
  
  # Annotate exposure/outcome
  mr_results[, `:=`(
    exposure = exposure_id,
    outcome = outcome_id,
    id.exposure = NULL,
    id.outcome = NULL
  )]
  setkey(mr_results, outcome)
  setkey(outcome_annotation, phenotype_id)
  mr_results <- mr_results[
    outcome_annotation[
      ,.(phenotype_id, outcome_long, subtype_long)], 
    nomatch = NULL
  ]
  
  # I need to know which variant set was used
  mr_results[, variant_set := variant_set_id]
  
  return(mr_results)
}
parallel::stopCluster(cl)

# SAVE ----

all_outcomes_mr_results <- rbindlist(
  all_outcomes_mr_results, 
  use.names = TRUE, 
  fill = TRUE
)
fwrite(
  all_outcomes_mr_results,
  glue("{result_save_path}/nmr_cvd_mr_results_{exposure_id}.tsv.gz")
)

# SESSION INFO ----

sessioninfo::session_info()
