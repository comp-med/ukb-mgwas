#!/bin/sh

# Running Regenie's step 1 for the EUR ancestry individuals

#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 30
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm_logs/step1.out
#SBATCH --chdir=/sc-projects/sc-proj-computational-medicine/people/Martijn/02_UKBB_NMR/03_GWAS_NMR/04_gwas_ancestry_specific/01_gwas_pop_csa
#SBATCH -J RegStep1_popCSA

# General setup for Regenie
regenie=/sc-projects/sc-proj-computational-medicine/programs/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl
tmpdir=tmp
prefix_out=output/preds/step1

# Genetic data that is needed
# For regenie step 1, that is the raw array data and the samples we want to include


# Do some additional filtering on the variants. Because we are using a subset of the variants, some of them have no variance anymore.

# 1. Create subset of the fam file for this subet of individuals
grep -f input/csa_samples_n9008_regenie.txt /sc-resources/ukb/data/bulk/genetic/array/step1_regenie_44448/pruned_genotypes/ukb22418_allChrs.pruned.fam > input/subset.fam

# 2. Filter on MAC for these individuals specifically
/sc-projects/sc-proj-computational-medicine/programs/plink/plink2_linux_amd_avx2_20240105/plink2 \
  --keep input/subset.fam \
  --bfile /sc-resources/ukb/data/bulk/genetic/array/step1_regenie_44448/pruned_genotypes/ukb22418_allChrs.pruned \
  --mac 100 \
  --write-snplist \
  --out output/snppass_mac100

# Pruned prefixes of the raw pruned data
pruned_genotypes=/sc-resources/ukb/data/bulk/genetic/array/step1_regenie_44448/pruned_genotypes/ukb22418_allChrs.pruned # Prefix of the pruned data
samples=input/csa_samples_n9008_regenie.txt                                                                             # Samples to keep

# Phenotypes / covariates
covariates=input/covariates.txt.gz
phenotypes=input/phenotypes.txt.gz

# Right locale settings needed for some reason
export LC_ALL=C \
export LANGUAGE=

${regenie} \
  --step 1 \
  --bed ${pruned_genotypes} \
  --keep ${samples} \
  --phenoFile ${phenotypes} \
  --covarFile ${covariates} \
  --threads 30 \
  --qt \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix ${tmpdir}/regenie_tmp_preds \
  --use-relative-path \
  --gz \
  --out ${prefix_out} \
  --extract output/snppass_mac100.snplist # The snp list we just created

echo '___________________'
echo 'Finished Step 1'

sbatch scripts/03_regenie_step2.sh
