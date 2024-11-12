rm(list = ls())
try(dev.off())

# Load required packages
library(ggplot2)
library(dplyr)
library(openxlsx)
library(magrittr)
library(reshape2)
library(data.table)
library(metafor)
library(patchwork)

setwd('')

mytheme <- readRDS('/home/mazo10/themes/theme.rds')


cs <- c(
  'HDL' = "#00468BFF",
  'LDL' = "#ED0000FF" ,
  'logTG' = "#42B540FF",
  'nonHDL' = "#0099B4FF",
  'TC' = "#925E9FFF"
)


# Get the LGC sentinel loci
lgc <- read.xlsx('input/41586_2021_4064_MOESM4_ESM.xlsx', sheet = 'ST.5 Multiancestry Results', rows = 5:2405)
lgc$markerName <- paste0('chr', lgc$chr, ':', lgc$`pos.(build.37)`) # We do this on absolute betas so no need for REF/ALT

lgc.traits <- unique(lgc$trait)

map <- openxlsx::read.xlsx('UK_biobank/nmr_telomer/qc/41597_2023_1949_MOESM2_ESM.xlsx') %>%
  filter(!is.na(UKB.Field.ID)) %>%
  filter(!Group %in% c('Amino acids', 'Ketone bodies', "Glycolysis related metabolites", "Fluid balance", "Inflammation", 'Fatty acids', 'Other lipids'))
nmr.traits <- map$Biomarker

# Get all the relevant loci from LGC
# Get own data that we need:
  # Absolute betas
  # SEs
  # LOG10Ps

# From this, generate:
  # ranking per LGC phenotype, per NMR phenotype of most strongly associated traits for these loci
  # Whether or not the phenotype is the 'intended' phenotype.

# ######
# LGC         NMR
# LDL-C       LDL_C - CHANGE TO CLINICAL_LDLC
# HDL-C       HDL_C
# TG          Total_TG
# Total C     Total_C
# nonHDL-C    non_HDL_C

# Conclusions:

# Among these loci, there are X instances where we do not find the matching / archetype trait among the top 10 most strong associations of the NMR platform

# Per trait, per locus
# Get the relevant stats for all NMR traits

lapply(lgc.traits, function(lgc.trait) {
  snps <- lgc %>% filter(trait == lgc.trait) %>% pull(rs_dbSNP150)
  snps <- snps[!is.na(snps)]

  # Get the stats
  ss <- bettermc::mclapply(nmr.traits, function(nmr.trait, snps) {

    fread(cmd = paste0("zcat 01_gwas/output/gwas_collated/gwas_", nmr.trait, "_allchr.txt.gz | grep -Ew 'CHROM|", paste0(snps, collapse = '|'), "'")) %>%
      mutate(pheno = nmr.trait)

  }, snps = snps, mc.cores = 25) %>%
    do.call(rbind, .) %>%
    mutate(lgctrait = lgc.trait)

  # Rankings of the most strongly associated
  ss %<>%
    mutate(abs_beta = abs(BETA)) %>%
    group_by(ID) %>%
    arrange(desc(LOG10P)) %>%
    mutate(ranking = 1:n()) %>%
    ungroup()

  return(ss)

}) %>% do.call(rbind, .) -> df

# Flag which NMR trait is the corresponding trait to the LGC trait
df$corresponding_trait <- F
df$corresponding_trait <- ifelse(df$lgctrait == 'HDL' & df$pheno == 'HDL_C', T, df$corresponding_trait)
df$corresponding_trait <- ifelse(df$lgctrait == 'LDL' & df$pheno == 'Clinical_LDL_C', T, df$corresponding_trait)
df$corresponding_trait <- ifelse(df$lgctrait == 'logTG' & df$pheno == 'Total_TG', T, df$corresponding_trait)
df$corresponding_trait <- ifelse(df$lgctrait == 'TC' & df$pheno == 'Total_C', T, df$corresponding_trait)
df$corresponding_trait <- ifelse(df$lgctrait == 'nonHDL' & df$pheno == 'non_HDL_C', T, df$corresponding_trait)


# Do some freaky merging to add the beta and SE of the corresponding trait as a column next to every row
df %>% filter(corresponding_trait) %>% mutate(B_corresponding = abs_beta, SE_corresponding = SE) %>% dplyr::select(c('ID', 'lgctrait', 'B_corresponding', 'SE_corresponding')) -> tmp
df <- merge(df, tmp, all.x=T, by = c('ID', 'lgctrait'))

# Meta-analysis on each row to see if the B and SE from the corresponding trait for that locus are in fact different
df$rownr <- 1:nrow(df)

bettermc::mclapply(df$rownr, function(ri) {

  betas = c(
    df %>% filter(rownr == ri) %>% pull(abs_beta),
    df %>% filter(rownr == ri) %>% pull(B_corresponding)
    )

  ses = c(
    df %>% filter(rownr == ri) %>% pull(SE),
    df %>% filter(rownr == ri) %>% pull(SE_corresponding)
  )

  ma <- rma(yi = betas, sei = ses, method = 'FE')

  # Get the correct P-value manually.
  data.frame(
    'ma_p' = ma$QEp,
    'rownr' = ri
  ) %>% return()

}, mc.cores = 15) %>% do.call(rbind, .) -> ma.res

df <- merge(df, ma.res, by = 'rownr')


# For the loci per trait, how often do we find the corresponding NMR trait in the top 10%?
# Top 10 of 205 traits corresponds to 20 -> e.g. how many times is the corresponding trait ranked <= 20?
df %>%
  mutate(top = ifelse(ranking <= 20, T, F)) %>%
  group_by(lgctrait) %>%
  filter(corresponding_trait) %>%
  summarise(incl_top = sum(top),
            not_incl_top = sum(top == F)) %>%
  mutate(totalloci = incl_top + not_incl_top) %>%
  mutate(frac_incl = incl_top / totalloci,
         frac_nonincl = not_incl_top / totalloci) -> stats

# Calculate an overall percentage of loci that have the corresponding stat in the top 10%?
df %>%
  mutate(top = ifelse(ranking <= 20, T, F)) %>%
  filter(corresponding_trait) %>%
  summarise(incl_top = sum(top),
            not_incl_top = sum(top == F)) %>%
  mutate(totalloci = incl_top + not_incl_top) %>%
  mutate(frac_incl = incl_top / totalloci,
         frac_nonincl = not_incl_top / totalloci)

# Get an idea of the 2166 loci, how many times do we find an NMR trait more sig associated compared to the archetype trait?
df$sig_diff <- ifelse(df$ma_p < (0.05 / 205), T, F)
df$high_eff <- ifelse(df$abs_beta > abs(df$B_corresponding), T, F)
df$sig <- ifelse(df$sig_diff & df$high_eff, T, F)


lapply(unique(df$lgctrait), function(tr) {

  tr.df <- df %>% filter(lgctrait == tr)

  bettermc::mclapply(unique(tr.df$ID) , function(x) {

    sub.df <- tr.df %>% filter(ID == x)

    return( data.frame(
      'locus' = x,
      'pos.corresp' = sub.df %>% filter(corresponding_trait) %>% pull(ranking),
      'any_sig_ma' = any(sub.df$sig)
    ))


  }, mc.cores = 20) %>% do.call(rbind, .) %>% as.data.frame() %>% mutate(pheno = tr)

}) %>% do.call(rbind, .) -> y

table(y$any_sig_ma) / nrow(y)



# Mutate the ranking so that the top rankings are on the right side
df$ranking_adj <- abs(df$ranking - 206)

pdf('graphics/distribution_ranks_lgctraits_correspondingNMRtrait.pdf', width = 2, height = 2)
ggplot(df %>% filter(corresponding_trait)) +
  geom_density(aes(x=ranking_adj, color = lgctrait, fill = lgctrait), alpha = .1) +
  mytheme +
  scale_color_manual(values = ggpubr::get_palette('lancet', 5)) +
  scale_fill_manual(values = ggpubr::get_palette('lancet', 5)) +
  labs(x = 'Ranking of NMR trait across loci',
       y = 'Density',
       fill = 'Global Lipids\nGenetics Consortium trait') +
  guides(colour = 'none') +
  geom_vline(xintercept = 185, lty = 2, color = 'lightgray') +
  theme(legend.position = 'none') +
  scale_x_continuous(breaks = c(6, 56, 106, 156, 206),
                     labels = function(x) {abs(x - 206)}) # Adjust the labels manually so that rank 1 is the highest, on the right side of the figure.
dev.off()

stats.melt <- melt(stats %>% select(-incl_top, -not_incl_top, -totalloci), id.vars = 'lgctrait')

pdf('graphics/piecharts_correspondingtrait_withintop10perc.pdf', width = 6, height = 2)
ggplot(stats.melt, aes(x = "", y = value, fill = lgctrait, color = lgctrait, alpha = variable)) +
  theme_void() +
  geom_col(width = 1) +
  coord_polar("y") +
  facet_wrap('lgctrait', nrow = 1) +
  geom_text(aes(label = paste0(round(value, 2) * 100, '%')), position = position_stack(vjust = 0.5), color = 'black') +
  scale_alpha_manual(values = c('frac_incl' = .4, 'frac_nonincl' = 0)) +
  scale_fill_manual(values = cs) +
  scale_color_manual(values = cs) +
  guides(alpha = 'none', fill = 'none', color = 'none')
dev.off()


# Results of the meta-analysis
# How often do we find traits that are sig. different from the corresponding trait?
# Correct the p-values for the nr of trait = 205, e.g. per locus and per phenotype
# alpha = 0.05 / 205
alpha = 0.05 / 205
df$sig <- ifelse(df$ma_p < alpha, T, F)
df$largerEff <- ifelse(df$abs_beta > df$B_corresponding, T, F)
df$sig2 <- ifelse(df$sig & df$largerEff, T, F)

# Per trait, per locus, how often do we find traits that are different?
df %>% group_by(lgctrait, ID) %>%
  summarise(sig = sum(sig2)) %>% arrange(desc(sig)) %>% head()


cs <- c(
  'gene' = 'forestgreen',
  'Triglycerides' = '#EE4C97',
  'Lipoprotein subclasses' = '#6F96AD',
  'Inflammation' = '#39806B',
  'Lipoprotein particle sizes' = '#7388AE',
  'Fatty acids' = '#917F59',
  'Relative lipoprotein lipid concentrations' = '#FCC292',
  'Amino acids' = '#BC3C29',
  'Fluid balance' = '#D58629',
  'Other lipids' = '#A1B0A3',
  'Cholesterol' = '#21689C',
  'Cholesteryl esters' = '#347693',
  'Free cholesterol' = '#868639',
  'Lipoprotein particle concentrations' = '#767AB0',
  'Phospholipids' = '#DDCC97',
  'Ketone bodies' = '#5E7A93',
  'Apolipoproteins' = '#6E5262',
  'Glycolysis related metabolites' = '#368549',
  'Total lipids' = '#F58794'
)
