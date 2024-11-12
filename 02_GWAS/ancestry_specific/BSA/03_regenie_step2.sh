#!/bin/sh

# Running Regenie's step 2

#SBATCH --partition=compute
#SBATCH --account=sc-users
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --array=1-23%23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 5
#SBATCH --output=slurm_logs/step2_%j.out
#SBATCH --chdir=/sc-projects/sc-proj-computational-medicine/people/Martijn/02_UKBB_NMR/03_GWAS_NMR/04_gwas_ancestry_specific/01_gwas_pop_csa
#SBATCH -J RegStep2_popCSA

# Chromosomes in parallel
chr=$(awk -v var=$SLURM_ARRAY_TASK_ID 'NR == var {print $1}' input/chromosomes.txt)

echo 'We are currently analysing chromosome:'
echo $chr

# General setup for Regenie
regenie=/sc-projects/sc-proj-computational-medicine/programs/regenie/regenie_v3.2.5.gz_x86_64_Centos7_mkl
tmpdir=tmp
regenie_step1_out=output/preds/step1_pred.list
out=output/gwas/chr${chr}

# Genetic data that is needed
bgendata=/sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c${chr}_b0_v3.bgen             # Imputed bgen files per chromosome
sample=/sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c${chr}_b0_v3.sample             # To check with Summaira

# Phenotypes / covariates
covariates=input/covariates.txt.gz
phenotypes=input/phenotypes.txt.gz
samples=input/csa_samples_n9008_regenie.txt                                                                             # Samples to keep

# Get the ancestry-specific variants to be included. Don't be too strict here, we can always filter out later.
# Filters:
  # HWE P ($8) > 1e-15
  # MAF ($14) > 0.1%

mkdir input/varqc_lists

# Necessary because the structure of chromosomes 1-22 and X is different
if [ $chr = "X" ]; then
    cat /sc-resources/ukb/data/bulk/genetic/imputed/variant_qc/CSA/ukb_imp_CSA_chr${chr}_snpstat.out | \
    grep -v '#' | \
    awk -F '\t' -v chr=${chr} ' { if ( (NR != 1 ) &&
                                     ( $8 > 1e-15 ) &&
                                     ( $17 > 0.001 ) )
                                     print $2 }' > input/varqc_lists/variants_chr${chr}.list

    echo "ChrX"
else

   cat /sc-resources/ukb/data/bulk/genetic/imputed/variant_qc/CSA/ukb_imp_CSA_chr${chr}_snpstat.out | \
   grep -v '#' | \
   awk -F '\t' -v chr=${chr} ' { if ( (NR != 1 ) &&
                                    ( $8 > 1e-15 ) &&
                                    ( $14 > 0.001 ) )
                                    print $2 }' > input/varqc_lists/variants_chr${chr}.list
fi

varqc_chr=input/varqc_lists/variants_chr${chr}.list

# Right locale settings needed for some reason
# Can also move this to bashrc, no need to do here
export LC_ALL=C \
export LANGUAGE=

${regenie} \
 --step 2 \
 --bgen  ${bgendata} \
 --keep ${samples} \
 --ref-first \
 --extract ${varqc_chr} \
 --sample ${sample} \
 --phenoFile ${phenotypes} \
 --covarFile ${covariates} \
 --threads 15 \
 --qt \
 --pred ${regenie_step1_out} \
 --bsize 1000 \
 --use-relative-path \
 --gz \
 --out ${out}
