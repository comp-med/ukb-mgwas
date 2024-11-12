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

# Script to annotate the sentinels variants for the EUR NMR GWAS


# Clumped GWAS results
lapply(fread('input/phenotypes.txt', header = F)$V1, function(x) {
  fread(paste0('output/sentinels/sentinels_', x, '.txt')) %>% mutate(pheno = x)
}) %>% do.call(rbind, .) %>%
  set_colnames(c('chrom', 'genpos', 'rsid', 'allele0', 'allele1', 'a1freq', 'info', 'n', 'test', 'beta', 'se', 'chisq', 'log10p', 'extra', 'begin', 'end', 'pheno')) %>%
  mutate(markerName = paste0('chr', chrom, ':', genpos, '_', pmin(allele0, allele1), '_', pmax(allele1, allele0))) %>%
  mutate(chrom = ifelse(chrom == 23, 'X', chrom)) -> df

# Regions indices
foo <- fread('input/inputtable.txt')
colnames(foo) <- c('chrom', 'genpos', 'pheno', 'rsid', 'region_name', 'proxy_job_index')
foo$chrom <- ifelse(foo$chrom == 23, 'X', foo$chrom)
df <- merge(df, foo, by = c('chrom', 'genpos', 'pheno', 'rsid'), all.x = F) # This removed 318 sentinels inside the MHC region


# Add the region information for good measure
regions <- fread('input/finemapping_regions.txt') %>%
  set_colnames(c('finemapping_idx', 'pheno', 'region_name', 'chrom', 'region_start', 'region_end')) %>%
  mutate(chrom = ifelse(chrom == 23, 'X', chrom))

df <- merge(df, regions, by = c('pheno', 'region_name', 'chrom'))

# Add the information by Nightingale
map <- openxlsx::read.xlsx('41597_2023_1949_MOESM2_ESM.xlsx')
df <- merge(df, map, by.x='pheno', by.y = 'Biomarker', all.x=T, all.y=F)


##########################################
####### LD proxy per sentinel
# Get the LD proxies per lead variant
lapply(paste0('output/proxies/', df$proxy_job_index, '/proxies.txt'), function(x) {
  fread(x) %>% filter(var != leadvar) %>% slice_max(order_by = R2, n = 1, with_ties = F)}) %>%
  do.call(rbind, .) %>%
  set_colnames(c('proxy', 'R2', 'leadvar', 'region_name', 'index')) -> ld.proxies

df <- merge(df, ld.proxies, by.x = c('proxy_job_index', 'region_name', 'rsid'), by.y = c('index', 'region_name', 'leadvar'), all.x=T)


##########################################
####### LD-clumping
LD_THR = 0.7

# Get proxies for all these sentinels
lapply(df$proxy_job_index, function(i, LD_THR) {
  fread(paste0('output/proxies/', i, '/proxies.txt')) %>% filter(R2 >= LD_THR)
}, LD_THR = LD_THR) %>% do.call(rbind, .) -> proxies

# Double check
summary(proxies$R2)

# All the sentinels should be there
# Double-check
stopifnot(all(df$rsid %in% proxies$leadvar))

# Order the in the correct way according to chromosome and genpos so that the first
# R2 group is actually R2group 1, 2, and so on from top down
df %>% arrange(chrom, genpos) %>% pull(proxy_job_index) -> ord

ld.sub     <- igraph::graph_from_data_frame(proxies %>% arrange(match(index, ord)) %>% select(leadvar, var, R2))
# Get the separate communities (= shared signals)
ld.sub     <- igraph::components(ld.sub)$membership
## convert to data frame
ld.sub     <- data.table(ID=names(ld.sub), R2.group=ld.sub)

# How many independent regions do we have?
summary(ld.sub$R2.group) # 2468 independently associated signals, missing MHC here
# Merge and arrange
df <- merge(df, ld.sub, by.x = 'rsid', by.y = 'ID') %>% arrange(chrom, genpos)


##########################################
###### 1. Add some of that Hg38 information
library(liftOver)
library(GenomicRanges)
print('Starting liftOver..')

# Unique list of variants to liftOver
df %>% pull(markerName) %>% unique() -> snp.lift

## subset
ukbb.snps  <- fread("QCed.UKBB.variants.Berlin.20220404.txt") %>%
  filter(MarkerName %in% snp.lift)

stopifnot(all(snp.lift %in% ukbb.snps$MarkerName)) # Do we have all the variants?

## edit chromosome
ukbb.snps$chr.hg39 <- ifelse(ukbb.snps$chromosome == 23, "X", ukbb.snps$chromosome) # Nah brah

## import chain
chain              <- import.chain("input/hg19ToHg38.over.chain")

## create GRanges object
grObject           <- GRanges(seqnames = paste0("chr", ukbb.snps$chr.hg39), ranges=IRanges(start = ukbb.snps$position, end = ukbb.snps$position, names=ukbb.snps$MarkerName))

## now do the mapping
tmp                <- as.data.frame(liftOver(grObject, chain))
## rename
names(tmp)         <- c("group", "MarkerName", "seqnames", "pos.hg38", "end", "width", "strand")

# Merge with df
df <- merge(df, tmp[, c('MarkerName', 'pos.hg38')], by.x='markerName', by.y='MarkerName', all.x=T)




##########################################
###### 2. GWAS Catalogue & openGWAS database to get previously annotated variants
print('Starting GWAS Catalogue')

# GWAS Catalogue
source('scripts/function_map_gwascatalogue.R')

# Read and prepare the GWAS catalogue
gwas.catalogue <- fread('input/gwas_catalog_v1.0.2-associations_e112_r2024-05-20.tsv', quote="")

## rename some columns
names(gwas.catalogue)[8] <- "TRAIT"

## prune GWAS catalogue data
gwas.catalogue           <- subset(gwas.catalogue, !is.na(`OR or BETA`) & is.finite(`OR or BETA`))
## generate risk allele and drop everything w/o this information
gwas.catalogue$riskA     <- sapply(gwas.catalogue$`STRONGEST SNP-RISK ALLELE`, function(x) strsplit(x,"-")[[1]][2])
## delete spaces at the end
gwas.catalogue$riskA     <- trimws(gwas.catalogue$riskA, which = "b")
## drop interaction entries
ii                       <- grep("[0-9]", gwas.catalogue$riskA)
gwas.catalogue           <- gwas.catalogue[-ii,]
## only genome-wide significant ones
gwas.catalogue           <- subset(gwas.catalogue, PVALUE_MLOG > 7.3)
## N = 420.745 entries, downloaded 14.05.2024

## create another entry to possible merge on (careful, genome build 38 mapping)
# gwas.catalogue[, snp.id := paste0(ifelse(CHR_ID == "X", 23, CHR_ID), ":", CHR_POS)]
# gwas.catalogue[, snp.id := CHR_ID, ":", CHR_POS]
gwas.catalogue %<>% mutate(snp.id = paste0(CHR_ID, ':', CHR_POS))

# To construct the new entry to match our lead variant and proxies on
ukbb.snps  <- fread("QCed.UKBB.variants.Berlin.20220404.txt")

# Get the proxies again to make sure these were not mutated somewhere above. We want all proxies
lapply(paste0('output/proxies/', df$proxy_job_index, '/proxies.txt'), function(x) { fread(x) }) %>%
  do.call(rbind, .) %>%
  set_colnames(c('proxy', 'R2', 'leadvar', 'region_name', 'index')) -> ld.proxies

# Skip for times sake
# Map our own data to the catalogue
# R2 = 0.8
gw.map.08 <- map.gwas.catalog(gwas.catalogue = gwas.catalogue, ukb.snps = ukb.snps, res.var = df, ld.proxies = ld.proxies, ld.thr = .8)
gw.map.08 <- gw.map.08 %>% distinct(leadvariant, .keep_all = T)

colnames(gw.map.08) <- paste0(colnames(gw.map.08), '_thrR2_0.8')
df <- merge(df, gw.map.08, by.x='rsid', by.y='leadvariant_thrR2_0.8', all.x=T)

# R2 = 0.1
gw.map.01 <- map.gwas.catalog(gwas.catalogue = gwas.catalogue, ukb.snps = ukb.snps, res.var = df, ld.proxies = ld.proxies, ld.thr = .1)
gw.map.01 <- gw.map.01 %>% distinct(leadvariant, .keep_all = T)

colnames(gw.map.01) <- paste0(colnames(gw.map.01), '_thrR2_0.1')
df <- merge(df, gw.map.01, by.x='rsid', by.y='leadvariant_thrR2_0.1', all.x=T)


##########################################
###### 5. VEP assignments

# Import VEP assignments that were prepared by Hannah
vep.snps   <- fread("VEP_QCed.UKBB.variants.Berlin.20222307_final.txt") %>%
  mutate(markerName = paste0('chr', CHROM, ':', POS, '_', pmin(REF, ALT), '_', pmax(REF, ALT)))

# Some primary variants are missing
table(df$markerName %in% vep.snps$markerName)

# Merge with the lead and proxy variants. Only the important proxies
imp.proxies <- ld.proxies %>% filter(R2 >= 0.8)

# Get the markerName in there
ukbb.snps  <- fread("QCed.UKBB.variants.Berlin.20220404.txt")
ukbb.snps %<>% dplyr::select(rsid, MarkerName)

imp.proxies <- merge(imp.proxies, ukbb.snps %>% set_colnames(c('rsid', 'leadvar_markerName')), by.x = 'leadvar', by.y = 'rsid', all.x=T, all.y=F)
imp.proxies <- merge(imp.proxies, ukbb.snps %>% set_colnames(c('rsid', 'proxy_markerName')), by.x = 'proxy', by.y = 'rsid', all.x=T, all.y=F)

# Still not al variants captured. It is what it is
table(unique(imp.proxies$leadvar_markerName) %in% vep.snps$markerName)
table(unique(imp.proxies$proxy_markerName) %in% vep.snps$markerName)


# merge the VEP
# for the main variant
imp.proxies <- merge(imp.proxies, vep.snps, by.x = 'leadvar_markerName', by.y = 'markerName', all.x=T, all.y=F)
# for the proxies
imp.proxies <- merge(imp.proxies, vep.snps, by.x = 'proxy_markerName', by.y = 'markerName', all.x=T, all.y=F, suffixes = c('.lead', '.proxy'))

## import consequence table
rank.cons  <- data.frame(consequence=c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
                                       "missense_variant", "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
                                       "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
                                       "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
                                       "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation", "regulatory_region_variant", "feature_truncation"),
                         rank=1:35)


imp.proxies$rank.proxy <- sapply(imp.proxies$Consequence.proxy, function(x){
  ## split multiple assignments
  x <- strsplit(x, "&")[[1]]
  # report most severe
  x <- subset(rank.cons, consequence %in% x)
  if(nrow(x) > 0){
    return(min(x$rank))
  }else{
    return(36)
  }
})

imp.proxies %<>% arrange(leadvar, rank.proxy) %>% group_by(leadvar) %>% mutate(ind = 1:n())
## Keep only the most severe proxy
sub <- imp.proxies %>% filter(ind == 1)
sub %<>% dplyr::select(c("leadvar_markerName","Consequence.lead", "IMPACT.lead", "SYMBOL.lead", "Gene.lead", "CADD_PHRED.lead",
                 "proxy", "R2", "POS.proxy", "ID.proxy", "Consequence.proxy", "IMPACT.proxy", "Gene.proxy",
                 "SYMBOL.proxy", "CADD_PHRED.proxy"))

# Merge with df
df <- merge(df, sub, by.x='markerName', by.y = 'leadvar_markerName', all.x=T)


##########################################
###### 6. Closest gene
require(biomaRt)

## get data on build 37
gene.ensembl                       <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)

## obtain a list of all protein coding genes
tmp.genes                          <- getBM(attributes = c('chromosome_name', 'start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'transcription_start_site'),
                                            filters = c('chromosome_name'),
                                            values = list(c(1:22, "X")),
                                            mart = gene.ensembl)

# Whats the stuff we are going to add?
df$closest.gene.body <- NULL; df$closest.gene.tss <- NULL; df$closest.genes <- NULL

df$closest.gene.body <- character()
df$closest.gene.tss <- character()
df$closest.genes <- character()

for (i in 1:nrow(df)) {
  tmp <- tmp.genes %>%
    filter(chromosome_name == df[[i, 'chrom']]) %>%
    filter(start_position >= df[[i, 'genpos']] - 1e6 & end_position <= df[[i, 'genpos']] + 1e6) %>%
    filter(gene_biotype %in% c("protein_coding", "processed_transcript"))

  # If results left, ..
  if(nrow(tmp) > 0) {
    # Calculate distance to gene bodies
    tmp %<>%
      rowwise() %>%
      mutate(dist_body = min(
        abs(start_position - df[[i, 'genpos']]),
        abs(end_position - df[[i, 'genpos']])
      )) %>%
      mutate(dist_tss = abs(transcription_start_site - df[[i, 'genpos']]))

    # Get top gene based on gene body distance and TSS distance
    top.genebody <- tmp %>% arrange(dist_body) %>% head(n=1)
    top.genetss <- tmp %>% arrange(dist_tss) %>% head(n=1)

    df$closest.gene.body[[i]] <- paste0(top.genebody$external_gene_name, ' (', top.genebody$dist_body, ') - ', top.genebody$ensembl_gene_id)
    df$closest.gene.tss[[i]] <- paste0(top.genetss$external_gene_name, ' (', top.genetss$dist_tss, ') - ', top.genetss$ensembl_gene_id)

    # Other closest genes within 0.5MB
    df$closest.genes[[i]] <- tmp %>% filter(dist_body < 5e5) %>% pull(external_gene_name) %>% unique() %>% paste0(., collapse = ' | ')
  }
}

print('All done, writing outputs')

# Export
# Complains about some cells that exceed the nr of permitted characters per cell
write.xlsx(df, file = 'output/annotated_loci.xlsx')
fwrite(df, file = 'output/annotated_loci.tsv', sep = '\t', quote=F, row.names = T, col.names = T)
