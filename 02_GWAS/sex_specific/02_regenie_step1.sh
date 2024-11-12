#!/bin/sh

# Regenie step 1 for males and females (EUR) separately

# General setup for Regenie
regenie=regenie_v3.2.5.gz_x86_64_Centos7_mkl

# Pruned prefixes of the raw pruned data
pruned_genotypes=ukb22418_allChrs.pruned # Prefix of the pruned data

tmpdir=tmpmales
prefix_out=output/preds/males/step1
samples=input/eur_males_n203300.txt                                                                                     # Samples to keep

# Right locale settings needed for some reason
export LC_ALL=C \
export LANGUAGE=



####### Males
# Phenotypes / covariates
covariates=input/covariates_males.txt
phenotypes=input/male_phenotypes.txt.gz

cmd="${regenie} \
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
  --out ${prefix_out}"

sbatch --partition=compute --account=sc-users --time=24:00:00 --nodes=1 --mem=20G \
  --ntasks=1 --cpus-per-task=30 --output=slurm_logs/males.step1.out \
  --chdir= \
  -J RegStep1_EURmales \
  --wrap "${cmd}"






####### Females
tmpdir=tmpfemales
prefix_out=output/preds/females/step1
samples=input/eur_females_n241027.txt                                                                                     # Samples to keep
# Phenotypes / covariates
covariates=input/covariates_females.txt
phenotypes=input/female_phenotypes.txt.gz

cmd="${regenie} \
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
  --out ${prefix_out}"

sbatch --partition=compute --account=sc-users --time=24:00:00 --nodes=1 --mem=20G \
  --ntasks=1 --cpus-per-task=30 --output=slurm_logs/females.step1.out \
  --chdir= \
  -J RegStep1_EURfemales \
  --wrap "${cmd}"
