#Try to close a graphics device if possible, and remove all from the env
rm(list = ls())
try(dev.off())

# Code from Aakash Nepal 

setwd('')

# extract the NMR biomarker data from the decoded data and give them short comprehensible column names
data <- data.table::fread("raw_data/data/20231213_category_220_221_222.tsv")
colnames(data) <- gsub("p(\\d+)_i(\\d+)", "f.\\1.\\2.0", colnames(data))

postqc <- ukbnmr::remove_technical_variation(data, version = 2)
saveRDS(postqc, file = 'output/postqc_nmr.rds')
