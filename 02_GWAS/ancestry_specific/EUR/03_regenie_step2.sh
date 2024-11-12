#!/bin/sh

# Running Regenie's step 2

#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --array=1-23%23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 35
#SBATCH --output=slurm_logs/step2_%j.out
#SBATCH --chdir=
#SBATCH -J RegStep2

# Chromosomes in parallel
chr=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR == var {print $1}' input/chromosomes.txt)

echo 'We are currently analysing chromosome:'
echo $chr

# General setup for Regenie
regenie=regenie_v3.2.5.gz_x86_64_Centos7_mkl
tmpdir=tmp
regenie_step1_out=output/preds/step1_pred.list
out=output/gwas/${chr}

# Genetic data that is needed
bgendata=imputed/bgen_files_44448/ukb22828_c${chr}_b0_v3.bgen             # Imputed bgen files per chromosome
sample=imputed/bgen_files_44448/ukb22828_c${chr}_b0_v3.sample             # To check with Summaira

# Phenotypes / covariates
covariates=input/covariates.txt.gz
phenotypes=input/phenotypes.txt.gz

# Get the variants to be included
# These are already QCed
mkdir input/varqc_lists
varqc=genotypes/variant_qc/output/ukb_imp_chr${chr}_qced.txt
cat ${varqc} | awk -v chr=${chr} '{if(NR != 1 && $3 == chr) print $2}' > input/varqc_lists/variants_chr${chr}.list
varqc_chr=input/varqc_lists/variants_chr${chr}.list

# Right locale settings needed for some reason
# Can also move this to bashrc, no need to do here
export LC_ALL=C \
export LANGUAGE=

${regenie} \
 --step 2 \
 --bgen  ${bgendata} \
 --ref-first \
 --extract ${varqc_chr} \
 --sample ${sample} \
 --phenoFile ${phenotypes} \
 --covarFile ${covariates} \
 --threads 35 \
 --qt \
 --pred ${regenie_step1_out} \
 --bsize 1000 \
 --use-relative-path \
 --gz \
 --out ${out}
