#Try to close a graphics device if possible, and remove all from the env
rm(list = ls())
try(dev.off())

# Load required packages
library(ggplot2)
library(dplyr)
library(openxlsx)
library(magrittr)
library(reshape2)
library(arrow)

# Working directory
setwd('')

# Get the imputed NMR data
load('working_data/working_data.Rdata')
nmr.phase12 <- data
sa.phase12 <- sample_annotation

load('working_data/working_data.Rdata')
nmr.phase3 <- data
sa.phase3 <- sample_annotation

# Remove the phase 1/2 variables that are not in phase3
nmr.phase12 %<>% select(-Corrected_Ala, -Glucose_Lactate)

stopifnot(all(colnames(nmr.phase12) == colnames(nmr.phase3)))

# Merge
nmr <- rbind(nmr.phase12, nmr.phase3)

cn <- intersect(colnames(sa.phase12), colnames(sa.phase3))
sa <- rbind(sa.phase12[, cn], sa.phase3[, cn])
sa$idx <- 1:nrow(sa)
# Get rid of trash
rm(data, nmr.phase12, nmr.phase3, sa.phase12, sa.phase3, sample_annotation, cn)


# Inverse rank normalise across both releases
# library(RNOmni)
# # Remove the eid column first, later re-insert
# apply(nmr %>% select(-eid), MARGIN = 2, RankNorm) -> nmr.irn
#
# nmr <- data.frame('eid' = nmr$eid, nmr.irn)


# Import the UKBB information
cl.select       <- c("f.eid", "f.21022.0.0", "f.31.0.0", "f.21001.0.0", 'f.6177.0.0', 'f.6177.1.0', 'f.6153.0.0', 'f.6153.1.0')
cl.names        <- c("f.eid", "age", "sex", "bmi", 'med_bl_male', 'med_fu_male', 'med_bl_female', 'med_fu_female')

## import data from the main release
ukb.dat         <- read_parquet("parquet_files/ukb45268.parquet",
                                col_select = cl.select)
## change names
names(ukb.dat)  <- cl.names

# Convert medication to characters
ukb.dat %<>% mutate(med_bl_male = as.character(med_bl_male), med_fu_male = as.character(med_fu_male),
                    med_bl_female = as.character(med_bl_female), med_fu_female = as.character(med_fu_female))

# Columns 6177 and 6153 present the most commonly prescribed medications
# Insulin, blood pressure, cholesterol: 6177, males
# Insulin, hormone replacement, contraception, blood pressure, cholesterol: 6153, females


# For both sexes, we fit a linear model in the form: lm(NMRfu vs NMRbl + age + sex + bmi + med1..medN)
######################
# Males

# Insulin is the main diabetic drug, but there are others: Metformin, etc, etc, etc that we ignore now.
# Correct for these manually and assess whether we have enough cases to estimate the effects
meduse <- data.table::fread('output/medication_usage_ukb.txt')

drugs <- c('Metformin')
bl.drugusers <- meduse %>% filter(drug.name %in% drugs) %>% filter(timepoint == 'Baseline') %>% filter(sex == 'Male') %>% pull(f.eid) %>% unique()
fu.drugusers <- meduse %>% filter(drug.name %in% drugs) %>% filter(timepoint == 'Followup') %>% filter(sex == 'Male') %>% pull(f.eid) %>% unique()

# Add this to the UKB dataframe
ukb.dat$med_bl_male <- ifelse(ukb.dat$med_bl_male == 'Insulin' | ukb.dat$f.eid %in% bl.drugusers, 'Diabetic medication', ukb.dat$med_bl_male)
ukb.dat$med_fu_male <- ifelse(ukb.dat$med_fu_male == 'Insulin' | ukb.dat$f.eid %in% fu.drugusers, 'Diabetic medication', ukb.dat$med_fu_male)

# Get the phenotypes of the people that did not take drugs a baseline.
phenotype.males <- ukb.dat %>%
  filter(sex =='Male') %>%
  mutate(med_bl_male = ifelse(med_bl_male == 'None of the above' | med_bl_male == 'Prefer not to answer' | med_bl_male == 'Do not know', 'No', med_bl_male)) %>%
  mutate(med_fu_male = ifelse(med_fu_male == 'None of the above' | med_fu_male == 'Prefer not to answer' | med_fu_male == 'Do not know', 'No', med_fu_male)) %>%
  filter(med_bl_male == 'No') # Filter for only the ones that took nothing at baseline

# We get only the individuals that have repeat measurements (e.g. appear exactly twice), and that did not take drugs at baseline as defined above
repeat.samples <- sa %>% filter(visit_index == 1) %>% pull(eid)

sa.males <- sa %>%
  filter(eid %in% phenotype.males$f.eid) %>% # Samples that did not take medication at bl
  filter(eid %in% repeat.samples) %>% # Samples that have repeat measurements
  group_by(eid) %>%
  add_count(name = "nr_visits") %>%
  filter(nr_visits == 2) %>%
  filter(!eid %in% c(1456300, 3048567, 4042425, 4449347))

table(sa.males$visit_index) # 6312 samples that have both timepoints and that did not take medication at baseline
phenotype.males %<>% filter(f.eid %in% sa.males$eid)
table(phenotype.males$med_fu_male)

# Out of which:
# Blood pressure medication     Cholesterol lowering medication   Diabetic medication           No
# 334                           776                               45                            5143

nmr.bl <- nmr[sa.males %>% filter(visit_index == 0) %>% pull(idx), ]
nmr.fu <- nmr[sa.males %>% filter(visit_index == 1) %>% pull(idx), ]

# Get everything in the right order
nmr.bl %<>% arrange(match(eid, phenotype.males$f.eid))
nmr.fu %<>% arrange(match(eid, phenotype.males$f.eid))
stopifnot(all(nmr.bl$eid == nmr.fu$eid) & all(nmr.fu$eid == phenotype.males$f.eid))

# Variables for linear model
age <- phenotype.males$age
bmi <- phenotype.males$bmi
med1 <- ifelse(phenotype.males$med_fu_male == 'Cholesterol lowering medication', 1, 0)
med2 <- ifelse(phenotype.males$med_fu_male == 'Blood pressure medication', 1, 0)
med3 <- ifelse(phenotype.males$med_fu_male == 'Diabetic medication', 1, 0)

res <- list()

for (metabolite in colnames(nmr.bl)) {
  if (metabolite == 'eid') {next}

  print(metabolite)
  met.bl <- nmr.bl[, metabolite]
  met.fu <- nmr.fu[, metabolite]

  res[[metabolite]] <- lm(met.fu ~ met.bl + age + bmi + med1 + med2 + med3) %>% broom::tidy() %>% mutate(met = metabolite)
}

do.call(rbind, res) %>%
  filter(term == 'med1' | term == 'med2' | term == 'med3') %>%
  mutate(term = recode_factor(term, 'med1' = 'Cholesterol-lowering medicine', 'med2' = 'Blood pressure medication', 'med3' = 'Diabetic medication')) -> male.res

write.table(x = male.res, file = 'output/all_medication/medication_metabolite_associations_followup_males.tsv', sep = '\t', quote = F, row.names = F, col.names = T)









#############################
# Females
# Insulin is the main diabetic drug, but there are others: Metformin, etc, etc, etc that we ignore now.
# Correct for these manually and assess whether we have enough cases to estimate the effects
meduse <- data.table::fread('output/medication_usage_ukb.txt')

drugs <- c('Metformin')
bl.drugusers <- meduse %>% filter(drug.name %in% drugs) %>% filter(timepoint == 'Baseline') %>% filter(sex == 'Female') %>% pull(f.eid) %>% unique()
fu.drugusers <- meduse %>% filter(drug.name %in% drugs) %>% filter(timepoint == 'Followup') %>% filter(sex == 'Female') %>% pull(f.eid) %>% unique()

# Add this to the UKB dataframe
ukb.dat$med_bl_female <- ifelse(ukb.dat$med_bl_female == 'Insulin' | ukb.dat$f.eid %in% bl.drugusers, 'Diabetic medication', ukb.dat$med_bl_female)
ukb.dat$med_fu_female <- ifelse(ukb.dat$med_fu_female == 'Insulin' | ukb.dat$f.eid %in% fu.drugusers, 'Diabetic medication', ukb.dat$med_fu_female)

# Get the phenotypes of the people that did not take drugs a baseline.
phenotype.females <- ukb.dat %>%
  filter(sex =='Female') %>%
  mutate(med_bl_female = ifelse(med_bl_female == 'None of the above' | med_bl_female == 'Prefer not to answer' | med_bl_female == 'Do not know', 'No', med_bl_female)) %>%
  mutate(med_fu_female = ifelse(med_fu_female == 'None of the above' | med_fu_female == 'Prefer not to answer' | med_fu_female == 'Do not know', 'No', med_fu_female)) %>%
  filter(med_bl_female == 'No') # Filter for only the ones that took nothing at baseline

# We get only the individuals that have repeat measurements (e.g. appear exactly twice), and that did not take drugs at baseline as defined above
repeat.samples <- sa %>% filter(visit_index == 1) %>% pull(eid)

sa.females <- sa %>%
  filter(eid %in% phenotype.females$f.eid) %>% # Samples that did not take medication at bl
  filter(eid %in% repeat.samples) %>% # Samples that have repeat measurements
  group_by(eid) %>%
  add_count(name = "nr_visits") %>%
  filter(nr_visits == 2) %>%
  filter(!eid %in% c(1199872)) # Remove one sample that has two timepoints for a repeat visit

table(sa.females$visit_index) # 6713 samples that have both timepoints and that did not take medication at baseline
phenotype.females %<>% filter(f.eid %in% sa.females$eid)
table(phenotype.females$med_fu_female) # Out of which 371 report cholesterol lowering medication at followup
# Blood pressure medication       Cholesterol lowering medication       Diabetic medication         Hormone replacement therapy        No       Oral contraceptive pill or minipill
# 257                             371                                   29                          148                                5881     27

nmr.bl <- nmr[sa.females %>% filter(visit_index == 0) %>% pull(idx), ]
nmr.fu <- nmr[sa.females %>% filter(visit_index == 1) %>% pull(idx), ]

# Get everything in the right order
nmr.bl %<>% arrange(match(eid, phenotype.females$f.eid))
nmr.fu %<>% arrange(match(eid, phenotype.females$f.eid))
stopifnot(all(nmr.bl$eid == nmr.fu$eid) & all(nmr.fu$eid == phenotype.females$f.eid))

# Variables for linear model
age <- phenotype.females$age
bmi <- phenotype.females$bmi
med1 <- ifelse(phenotype.females$med_fu_female == 'Cholesterol lowering medication', 1, 0)
med2 <- ifelse(phenotype.females$med_fu_female == 'Oral contraceptive pill or minipill', 1, 0)
med3 <- ifelse(phenotype.females$med_fu_female == 'Hormone replacement therapy', 1, 0)
med4 <- ifelse(phenotype.females$med_fu_female == 'Blood pressure medication', 1, 0)
med5 <- ifelse(phenotype.females$med_fu_female == 'Diabetic medication', 1, 0)

res <- list()

for (metabolite in colnames(nmr.bl)) {
  if (metabolite == 'eid') {next}

  print(metabolite)
  met.bl <- nmr.bl[, metabolite]
  met.fu <- nmr.fu[, metabolite]

  res[[metabolite]] <- lm(met.fu ~ met.bl + age + bmi + med1 + med2 + med3 + med4 + med5) %>% broom::tidy() %>% mutate(met = metabolite)
}

do.call(rbind, res) %>%
  filter(term == 'med1' | term == 'med2' | term == 'med3' | term == 'med4' | term == 'med5') %>%
  mutate(term = recode_factor(term,
                              'med1' = 'Cholesterol-lowering medicine',
                              'med2' = 'Oral contraceptive pill or minipill',
                              'med3' = 'Hormone replacement therapy',
                              'med4' = 'Blood pressure medication',
                              'med5' = 'Diabetic medication')) -> female.res

write.table(x = female.res, file = 'output/all_medication/medication_metabolite_associations_followup_females.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
