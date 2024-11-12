rm(list = ls())
try(dev.off())

# Load required packages
library(ggplot2)
library(dplyr)
library(openxlsx)
library(magrittr)
library(reshape2)
library(data.table)

setwd('')

mytheme <- readRDS('/home/mazo10/themes/theme.rds')

phenotypes <- fread('../01_gwas/input/phenotypes_mapped.txt', header=F)$V1

lapply(phenotypes, function(phenotype) {
  dir.create(paste0('output/', phenotype))

  # Get top varirs
  topvars <- fread('Results.NMR.lead.credible.set.variant.annotation.20240719.txt') %>% 
    filter(pheno.locus == phenotype)

  # Write files per chr on what snps to take
  lapply(unique(topvars$chrom), function(chr) {

    topvars %>% filter(chrom == chr) %>% pull(id) %>% cat(., sep = '\n', file = paste0('output/', phenotype, '/', chr, '_snps.txt'))

  })

  # Get the dosages
  lapply(unique(topvars$chrom), function(chr) {

    cmd = paste0('sh scripts/02_get_variants.sh ', phenotype, ' ', chr)
    system(cmd)

  })

  # Read the dosages in
  lapply(unique(topvars$chrom), function(chr) {

    fread(paste0('output/', phenotype, '/snps_', chr, '_dosage.txt'))

  }) %>% rbindlist(., fill=T) -> dosage.mat

  dosage.mat %<>% tibble::column_to_rownames('SNPID') %>%
    select(-chromosome, -rsid, -position, -alleleA, -alleleB) %>%
    t()

  # Get the phenotypes
  fread('../01_gwas/input/phenotypes.txt.gz') %>% tibble::column_to_rownames('FID') %>% select(phenotype) %>% as.matrix() -> phenotypes


  # Same samples, make sure
  samples <- intersect(rownames(dosage.mat), rownames(phenotypes))

  dosage.mat[samples, ] -> dosage.mat
  phenotypes[samples, ] %>% as.matrix() %>% set_colnames(c('nmrpheno'))-> phenotypes

  stopifnot(all(rownames(phenotype) == rownames(dosage.mat)))

  testdf <- cbind(dosage.mat, phenotypes)

  # Linear model to model all these variants on metabolite levels.
  # Report adjusted R2 per metabolite from combined model

  # dosage.mat$nmrpheno <- phenotypes

  lm(data = as.data.frame(testdf), formula = nmrpheno ~ .) -> x

  r2 <- summary(x)$adj.r.squared

  out <- data.frame(
    pheno = phenotype,
    adjR2 = r2
  )

  return(out)
}) %>% do.call(rbind, .) -> tmp

fwrite(tmp, 'output/r2_perphenotype_collated.tsv', sep = '\t')
