#Try to close a graphics device if possible, and remove all from the env
rm(list = ls())
try(dev.off())

# Load required packages
library(dplyr)
library(magrittr)
library(mice)

setwd('')

# Read the data
nmr <- read.table('output/postqc_cleaned/data.tsv', sep='\t', header=T)
fa <- read.table('output/postqc_cleaned/feature_annotation.tsv', sep='\t', header=T)
sa <- read.table('output/postqc_cleaned/sample_annotation.tsv', sep='\t', header=T)

# Get rid of the QC FAIL samples before imputation
sa$idx <- 1:nrow(sa)
sa %<>% filter(exclude_qc == FALSE)
nmr <- nmr[sa$idx, ]

stopifnot(all(sa$eid == nmr$eid))


# Impute the nmr measurements
# In total, we are imputing: 0.16% of data
( sum(is.na(nmr)) / prod(dim(nmr)) ) * 100

to.imp <- nmr %>% select(-eid, -visit_index)

pred <- quickpred(to.imp, mincor=.3)

imp <- mice(to.imp,
            predictorMatrix = pred,
            method = 'pmm',
            m = 1,
            remove.collinear=F)

write.table(imp$loggedEvents, file = 'slurm_logs/imputation_loggedevents.txt', sep='\t', quote=F)

imp <- complete(imp)
imp <- cbind(data.frame(eid=sa$eid), imp)

# Export the tables
write.table(imp, file = 'output/imputed/data_imputed.tsv', sep = '\t', eol = '\n', quote = F, row.names = F, col.names = T)
write.table(fa, file = 'output/imputed/feature_annotation.tsv', sep = '\t', eol = '\n', quote = F, row.names = F, col.names = T)
write.table(sa, file = 'output/imputed/sample_annotation.tsv', sep = '\t', eol = '\n', quote = F, row.names = F, col.names = T)
