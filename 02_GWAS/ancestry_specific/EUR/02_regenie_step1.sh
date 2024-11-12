#!/bin/sh

# Running Regenie's step 1 for the EUR ancestry individuals

#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 30
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm_logs/step1.out
#SBATCH --chdir=
#SBATCH -J RegStep1

# General setup for Regenie
regenie=regenie_v3.2.5.gz_x86_64_Centos7_mkl
tmpdir=tmp
prefix_out=output/preds/step1

# Genetic data that is needed
# For regenie step 1, that is the raw array data and the samples we want to include


# Pruned prefixes of the raw pruned data
pruned_genotypes=step1_regenie_44448/pruned_genotypes/ukb22418_allChrs.pruned # Prefix of the pruned data
samples=input/european_samples_n444284_regenie.txt                                                                      # Samples to keep

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
  --out ${prefix_out}
