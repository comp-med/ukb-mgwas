#-------------------------------------------------------------------------#
# Metabolic QTL Data Processing Pipeline
# Main steps:
# 1. Load and process variants (both lead and proxies)
# 2. Map coordinates between genome builds
# 3. Analyze gene context and LD regions
# 4. Add functional annotations (eQTLs, VEP, pathways)
# 5. Create pathway-specific datasets
#-------------------------------------------------------------------------#

## Required libraries
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(data.table)
library(glue)
library(tidyverse)
library(magrittr)

## Source helper functions
source("settings_classifier.R")
source("functions_classifier.R")

#-------------------------------------------------------------------------#
# Step 1: Load and Process Variants
#-------------------------------------------------------------------------#

## Load lead variants from fine-mapping
finemapped_results_NMR_top_var <- fread("lead_variants_finemapped.txt.gz", header = TRUE)
combined_lead_data <- fread("collated_finemapped_results.txt.gz", header = TRUE)
combined_lead_data <- combined_lead_data[, .(chrom, genpos, id, allele0, allele1)]
combined_lead_data <- unique(combined_lead_data)

## Load and process proxy variants (r2 > 0.6)
proxy_files <- list.files("proxies_dir", full.names = TRUE)
proxy_data_list <- lapply(proxy_files, fread)
proxy_combined_data <- rbindlist(proxy_data_list, use.names = TRUE, fill = TRUE)
proxy_combined_data <- proxy_combined_data[r2 > 0.6]

#-------------------------------------------------------------------------#
# Step 2: Coordinate Mapping (hg19 to hg38)
#-------------------------------------------------------------------------#

## Map coordinates using liftOver and biomart as backup
chain <- import.chain("hg19ToHg38.over.chain")
grObject_leads <- GRanges(
  seqnames = paste0("chr", combined_lead_data$chr.hg38),
  ranges = IRanges(
    start = combined_lead_data$pos,
    end = combined_lead_data$pos
  )
)
mapped_coordinates <- liftOver(grObject_leads, chain)

#-------------------------------------------------------------------------#
# Step 3: Gene Context Analysis
#-------------------------------------------------------------------------#

## Load gene annotations
granges_38 <- import("genome.gtf", format = "gtf")

## Find overlapping genes within defined windows
local_genes <- find_overlapping_genes(lead_ranges, granges_38, interval_kb)
local_genes <- gene_variant_distance_finder(lead_ranges, local_genes)
nearest_genes<- optimized_nearest_genes_selector(
  local_genes, local_genes_annotation, biotype_of_interest,
  NUMBER_NEAREST
)
## Analyze LD regions
LD_region_ranges <- LD_region_range_finder(lead_data, proxy_data)
LD_region_genes <- find_overlapping_genes(LD_region_ranges, granges_38, LD_region_overhang_kb)

classifier_dataset  <- merge(
  nearest_genes,
  LD_region_genes,
  by.x = c("ensembl_id", "LEAD_rsID"),
  by.y = c("ensembl_id", "LEAD_rsID"),
  all = TRUE
)


#-------------------------------------------------------------------------#
# Step 4: Add Functional Annotations
#-------------------------------------------------------------------------#

#-------------------------------------------------------------------------#
## 4.1 Process eQTL Data from eQTL Catalogue
# - Filter for genome-wide significant associations (p < 5e-8)
# - Create tissue-specific matrices for each variant-gene pair
#-------------------------------------------------------------------------#

## Read and process eQTL data
eqtl_catalogue_files <- list.files(pattern = "\\.tsv\\.gz$")
eqtl_data <- lapply(eqtl_catalogue_files, function(filepath) {
  # Read essential columns
  data <- fread(filepath, select = c("variant", "rsid", "gene_id", "pvalue", "pip"))
  # Filter for significance
  data <- data[pvalue <= 5e-8]
  return(data)
})

## Create tissue-specific PIP matrix
eqtl_catalogue <- rbindlist(eqtl_data)
# Take highest PIP for each tissue-variant-gene combination
eqtl_matrix <- eqtl_catalogue[
  , .SD[which.max(pip)],
  by = .(tissue_label, variant_id, gene_id)
]
# Pivot to create tissue-specific columns
eqtl_matrix_wide <- dcast(
  eqtl_matrix,
  variant_id + gene_id ~ tissue_label,
  value.var = "pip"
)

## Combine all annotations
classifier_dataset <- merge(
  classifier_dataset,
  eqtl_matrix_wide,
  by.x = c("ensembl_id", "variant_id"),
  by.y = c("gene_id", "variant_id"),
  all.x = TRUE
)

#-------------------------------------------------------------------------#
## 4.2 VEP Annotation Integration
# - Focus on HIGH and MODERATE impact variants
# - Process both lead and proxy variants
# - Integrate most severe consequence per gene
#-------------------------------------------------------------------------#

## Process VEP annotations
VEP_annotation <- VEP_annotation[
  IMPACT %in% c("HIGH", "MODERATE")
]

## Find most severe consequence for each variant-gene pair
severity_ranking <- c(
  "transcript_ablation" = 1,
  "splice_acceptor_variant" = 2,
  # ... other consequences in order of severity
  "regulatory_region_variant" = 33
)

## Get most severe consequence per variant
VEP_processed <- VEP_annotation[
  , .(consequence = consequence[which.min(severity_ranking[consequence])]),
  by = .(variant_id, gene_id)
]

## Add VEP annotations
classifier_dataset <- merge(
  classifier_dataset,
  VEP_processed,
  by.x = c("ensembl_id", "variant_id"),
  by.y = c("gene_id", "variant_id"),
  all.x = TRUE
)

#-------------------------------------------------------------------------#
## 4.3 Add approved drug targets, OMIM annotations, MGI and ORPHANET annotations from PROGEM (metabolic genes)
#-------------------------------------------------------------------------#
approved_drug_targets <- fread("approved_drug_targets.tsv", sep = "\t")

classifier_dataset$approved_drug_target <-
  classifier_dataset$ensembl_id %in% approved_drug_targets_$id

omim_genes_catalogue <- fread("omim/mim2gene.txt")
classifier_dataset$omim <-
  classifier_dataset$ensembl_id %in% omim_genes_catalogue$`Ensembl Gene ID (Ensembl)`
classifier_dataset$omim <- as.numeric(classifier_dataset$omim)

mgi_orphanet_annotations <- fread("mgi_orphanet_progem", sep = "\t")
classifier_dataset$mgi_annotation <-
  classifier_dataset$ensembl_id %in% mgi_orphanet_annotations$mgi_gene
classifier_dataset$orphanet_annotation <-
  classifier_dataset$ensembl_id %in% mgi_orphanet_annotations$orphanet_gene

#-------------------------------------------------------------------------#
## 4.4 Pathway-Specific Dataset Creation
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
### Prepare Reactome Annotation ---------------------------------------------
#-------------------------------------------------------------------------#
require(ReactomeContentService4R)
require(parallel)
## get pathways
reactome.pathway <- getSchemaClass(class = "Pathway", species = "human", all = TRUE)
## drop disease pathways
reactome.pathway <- as.data.table(reactome.pathway)
reactome.pathway <- reactome.pathway[isInDisease == F]
## drop top level pathways
reactome.pathway <- reactome.pathway[schemaClass == "Pathway"]
## get associated genes
reactome.pathway <- mclapply(1:nrow(reactome.pathway), function(x) {
  ## obtain relevant genes
  tmp <- event2Ids(reactome.pathway$stId[x])
  ## return information
  return(data.table(reactome.pathway[x, ], gene.symbol = tmp$geneSymbol))
}, mc.cores = 10)
## combine
reactome.pathway <- rbindlist(reactome.pathway, fill = T)
reactome.pathway[, name := NULL]
## make unique
reactome.pathway <- unique(reactome.pathway)
## exclude some very large, possibly unspecific pathways
tail(sort(table(reactome.pathway$displayName)))

#-------------------------------------------------------------------------#
### Prepare KEGG Annotation -------------------------------------------------
#-------------------------------------------------------------------------#
library(limma)
### We get entrez ids and their pathways.
kegg.pathway <- getGeneKEGGLinks(species = "hsa")
### This is to get the gene symbols using entrez ids
kegg.pathway$Symbol <- mapIds(org.Hs.eg.db, kegg.pathway$GeneID, column = "SYMBOL", keytype = "ENTREZID")
## pathway names
pathway_names <- getKEGGPathwayNames(species = "hsa")
### obtain all human KEGG pathways
kegg.pathway <- merge(kegg.pathway, pathway_names, by = "PathwayID")
kegg.pathway <- as.data.table(kegg.pathway)
### drop some very large, but unspecific pathways
kegg.pathway <- kegg.pathway[!(Description %in% c(
  "Metabolic pathways - Homo sapiens (human)", "Pathways in cancer - Homo sapiens (human)",
  "Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)", "MicroRNAs in cancer - Homo sapiens (human)",
  "Coronavirus disease - COVID-19 - Homo sapiens (human)", "Chemical carcinogenesis - reactive oxygen species - Homo sapiens (human)",
  "Proteoglycans in cancer - Homo sapiens (human)", "Chemical carcinogenesis - receptor activation - Homo sapiens (human)",
  "Transcriptional misregulation in cancer - Homo sapiens (human)"
))]

#-------------------------------------------------------------------------#
### 4.4.1 Cholesterol Metabolism Dataset
#-------------------------------------------------------------------------#

## Define relevant pathways
cholesterol_kegg_pathways <- c(
  "hsa04979", # Cholesterol metabolism
  "hsa00100", # Steroid biosynthesis
  "hsa00984", # Steroid degradation
  "hsa00120", # Primary bile acid biosynthesis
  "hsa04975", # Fat digestion and absorption
  "hsa04976", # Bile secretion
  "hsa00900" # Terpenoid backbone biosynthesis
)

cholesterol_reactome_pathways <- c(
  "R-HSA-191273", # Cholesterol biosynthesis
  "R-HSA-6807062", # Cholesterol biosynthesis via lathosterol
  "R-HSA-6807047", # Cholesterol biosynthesis via desmosterol
  "R-HSA-1655829", # Regulation of cholesterol biosynthesis by SREBP (SREBF)
  "R-HSA-8963678", # Intestinal lipid absorption
  "R-HSA-193807", # Synthesis of bile acids and bile salts via 27-hydroxycholesterol
  "R-HSA-193775", # Synthesis of bile acids and bile salts via 24-hydroxycholesterol
  "R-HSA-196071", # Metabolism of steroid hormones
  "R-HSA-194002", # Glucocorticoid biosynthesis
  "R-HSA-193144", # Estrogen biosynthesis
  "R-HSA-193993", # Mineralocorticoid biosynthesis
  "R-HSA-192105", # Synthesis of bile acids and bile salts
  "R-HSA-8964058", # HDL remodeling
  "R-HSA-2426168", # Activation of gene expression by SREBF (SREBP)
  "R-HSA-211976", # Endogenous sterols
  "R-HSA-8963896", # HDL assembly
  "R-HSA-194068", # Bile acid and bile salt metabolism
  "R-HSA-1369062", # ABC transporters in lipid homeostasis
  "R-HSA-2173782", # Binding and Uptake of Ligands by Scavenger Receptors
  "R-HSA-8964041", # LDL remodeling
  "R-HSA-192456", # Digestion of dietary lipid
  "R-HSA-8963888", # Chylomicron assembly
  "R-HSA-8964038", # LDL clearance
  "R-HSA-8979227", # Triglyceride metabolism
  "R-HSA-8866423", # VLDL assembly
  "R-HSA-8963899", # Plasma lipoprotein remodeling
  "R-HSA-163560", # Triglyceride catabolism
  "R-HSA-174824", # Plasma lipoprotein assembly, remodeling, and clearance
  "R-HSA-8964043" # Plasma lipoprotein clearance
)

## Create cholesterol-specific dataset
cholesterol_dataset <- function(classifier_dataset) {
  # Get cholesterol-related variants
  cholesterol_variants <- annotation_phenotypes[
    Group %in% c("Cholesterol", "Cholesteryl esters", "Free cholesterol")
  ]

  # Subset main dataset
  subset_classifier <- classifier_dataset[
    LEAD_rsID %in% cholesterol_variants$rsID
  ]

  # Add pathway annotations
  subset_classifier[
    , kegg_pathway_specific := as.integer(
      hgnc_symbol %in% cholesterol_pathways_kegg$Symbol
    )
  ][
    , reactome_pathway_specific := as.integer(
      hgnc_symbol %in% reactome_pathways_cholesterol$gene.symbol
    )
  ]

  return(subset_classifier)
}

#-------------------------------------------------------------------------#
### 4.4.2 Amino Acid Metabolism Dataset
#-------------------------------------------------------------------------#

## Define relevant pathways
amino_acid_kegg_pathways <- c(
  "hsa01230", # Biosynthesis of amino acids
  "hsa00260", # Glycine, serine and threonine metabolism
  "hsa00290", # Valine, leucine and isoleucine biosynthesis
  "hsa00400", # Phenylalanine, tyrosine and tryptophan biosynthesis - hsa00400
  "hsa00360", # Phenylalanine metabolism
  "hsa00340", # Histidine metabolism - hsa00340
  "hsa00280", # Valine, leucine and isoleucine degradation - hsa00280
  "hsa00350", # Tyrosine metabolism - hsa00350
  "hsa00250", # Alanine, aspartate and glutamate metabolism - hsa00250
  "hsa00020" # Citrate cycle (TCA cycle) - hsa00020
)

amino_acid_reactome_pathways <- c(
  "R-HSA-71291", # Metabolism of amino acids and derivatives -
  "R-HSA-352230", # Amino acid transport across the plasma membrane
  "R-HSA-70635", # Urea cycle
  "R-HSA-70895", # Branched-chain amino acid catabolism
  "R-HSA-8964540", # Alanine metabolism
  "R-HSA-8964539", # Glutamate and glutamine metabolism
  "R-HSA-8963691", # Phenylalanine and tyrosine metabolism
  "R-HSA-6783984", # Glycine degradation
  "R-HSA-70268", # Pyruvate metabolism (because alanine)
  "R-HSA-8963684", # Tyrosine catabolism
  "R-HSA-8964208", # Phenylalanine metabolism
  "R-HSA-70921" #  Histidine catabolism
)

## Create amino acid dataset using similar approach to cholesterol dataset

#-------------------------------------------------------------------------#
### 4.4.3 Lipid Metabolism Dataset
#-------------------------------------------------------------------------#

## Define relevant pathways
lipid_kegg_pathways <- c(
  "hsa04975", # Fat digestion and absorption
  "hsa02010", # ABC transporters
  "hsa00071", # Fatty acid degradation
  "hsa00061", # Fatty acid biosynthesis
  "hsa00062", # Fatty acid elongation
  "hsa01040", # Biosynthesis of unsaturated fatty acids
  "hsa01212", # Fatty acid metabolism
  "hsa04923", # Regulation of lipolysis in adipocytes
  "hsa00600", # Sphingolipid metabolism
  "hsa00561" #  Glycerolipid metabolism
)

lipid_reactome_pathways <- c(
  "R-HSA-556833", # Metabolism of lipids
  "R-HSA-8963898", # Plasma lipoprotein assembly
  "R-HSA-8978868", # Fatty acid metabolism
  "R-HSA-1483257", # Phospholipid metabolism
  "R-HSA-9840310", # Glycosphingolipid catabolism
  "R-HSA-8964041", # LDL remodeling
  "R-HSA-8964572", # Lipid particle organization
  "R-HSA-390918", # Peroxisomal lipid metabolism
  "R-HSA-8963678", # Intestinal lipid absorption
  "R-HSA-1369062", # ABC transporters in lipid homeostasis
  "R-HSA-1483206", # Glycerophospholipid biosynthesis
  "R-HSA-174824", # Plasma lipoprotein assembly, remodeling, and clearance
  "R-HSA-8866423", # VLDL assembly
  "R-HSA-8963896", # HDL assembly
  "R-HSA-8963899", # Plasma lipoprotein remodeling
  "R-HSA-163560", # Triglyceride catabolism
  "R-HSA-8979227" # Triglyceride metabolism
)

## Create lipid dataset using similar approach to cholesterol dataset


## Write final datasets
fwrite(classifier_dataset, "output/classifier_dataset.tsv")
fwrite(cholesterol_dataset, "output/classifier_dataset_cholesterol.tsv")
fwrite(amino_acid_dataset, "output/classifier_dataset_amino_acids.tsv")
fwrite(lipid_dataset, "output/classifier_dataset_lipids.tsv")
