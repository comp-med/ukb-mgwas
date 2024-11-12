#Try to close a graphics device if possible, and remove all from the env
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

# Script to collate the regions that were defined earlier and prepare the input for Susie.

## collate results files from regional clumping
list.files('output/sentinels/', full.names = T) -> files
# files <- files[!grepl('unadjusted', files)]

res <- lapply(files, function(f) {

  fread(f) %>%
    set_colnames(c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "region_start", "region_end")) %>%
    mutate(phenotype = gsub(pattern = 'output/sentinels//sentinels_', replacement = '', f, fixed = T)) %>%
    mutate(phenotype = gsub(pattern = '.txt', '', phenotype, fixed=T)) %>%
    arrange(CHROM, region_start)

  }) %>% do.call(rbind, .)

# Narrow the regions
res$region_start_cap <- res$region_start + 2e5
res$region_end_cap   <- res$region_end - 2e5

res %<>% arrange(CHROM, region_start_cap)


# Taken from: https://stackoverflow.com/questions/28938147/how-to-flatten-merge-overlapping-time-periods
collapsed.regions <- res %>%
  arrange(CHROM, region_start_cap) %>%
  group_by(CHROM) %>%
  mutate(indx = c(0,
                  cumsum(as.numeric(lead(region_start_cap)) >
                           cummax(as.numeric(region_end_cap)))[-n()])) %>%
  group_by(CHROM, indx) %>%
  summarise(region_start_cap = data.table::first(region_start_cap), region_end_cap = last(region_end_cap)) %>%
  as.data.table

setkey(collapsed.regions, CHROM, region_start_cap, region_end_cap)

# The prefix 'i' represent the former regional boundaries before merging with the newly defined regions in 'collapsed.regions'
res <- foverlaps(res, collapsed.regions, by.x = c('CHROM', "region_start_cap", 'region_end_cap'))

# New index to denote the regions. Format: CHROM_locusnr
res$region <- paste0(res$CHROM, '_', res$indx)

length(unique(res$region))

miss <- ifelse(res$GENPOS < res$region_end_cap & res$GENPOS > res$region_start_cap, T, F)
table(miss)

# We are missing signals, expand these regions only
res$miss <- miss
res %>% filter(miss==F) %>% pull(region) -> adj
res$adj_region <- ifelse(res$region %in% adj, T, F)

res %<>%
  mutate(region_start_cap = ifelse(!adj_region, region_start_cap, region_start_cap - 10e5)) %>%
  mutate(region_end_cap = ifelse(!adj_region, region_end_cap, region_end_cap + 10e5))

# Edge cases
res$region_start_cap <- ifelse(res$region_start_cap < 0, 0, res$region_start_cap)

miss <- ifelse(res$GENPOS < res$region_end_cap & res$GENPOS > res$region_start_cap, T, F)
table(miss) # We cover all signals now

res$width <- res$region_end_cap - res$region_start_cap

# Identify the MHC region and get that out for now.
# chr6:25.5-34.0Mb
res %<>% mutate(mhc = ifelse(CHROM == 6 & GENPOS > 25.5e6 & GENPOS < 34.0e6, T, F))
table(res$mhc) # Sentinel variants in the MHC region
res %>% filter(mhc) %>% pull(region) %>% unique() # Region 6_10 is the MHC one, remove
res %<>% filter(region != '6_10')

# Clean up
res %<>%
  select(-i.region_end_cap, -i.region_start_cap, -EXTRA, -region_start, -region_end)

# Input for Susie
res %>%
  distinct(phenotype, region, CHROM, region_start_cap, region_end_cap, width) %>%
  mutate(idx = 1:n()) %>%
  select(idx, phenotype, region, CHROM, region_start_cap, region_end_cap) -> regions

write.table(regions, file = 'input/finemapping_regions.txt', sep="\t", col.names = F, row.names = F, quote = F)
