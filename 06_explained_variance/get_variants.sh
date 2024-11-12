#!/bin/sh

## export location of files
export dir=bgen_files_44448/

phenotype=$1
chr=$2


if [ ${chr} -eq 23 ]; then

  ## create subset bgen file get only SNPs in the data set
  programs/bgen/build/apps/bgenix \
    -g ${dir}/ukb22828_cX_b0_v3.bgen \
    -incl-rsids output/${phenotype}/${chr}_snps.txt > output/${phenotype}/flt_dosages.bgen

  programs/qctool/build/release/qctool_v2.0.7 \
    -threads 5 \
    -g output/${phenotype}/flt_dosages.bgen \
    -s ${dir}/ukb22828_cX_b0_v3.sample \
    -incl-samples input/european_samples_n444284_regenie.txt \
    -og - \
    -ofiletype dosage > output/${phenotype}/snps_${chr}_dosage.txt

else

  ## create subset bgen file get only SNPs in the data set
  programs/bgen/build/apps/bgenix \
    -g ${dir}/ukb22828_c${chr}_b0_v3.bgen \
    -incl-rsids output/${phenotype}/${chr}_snps.txt > output/${phenotype}/flt_dosages.bgen

programs/qctool/build/release/qctool_v2.0.7 \
    -threads 5 \
    -g output/${phenotype}/flt_dosages.bgen \
    -s ${dir}/ukb22828_c${chr}_b0_v3.sample \
    -incl-samples input/european_samples_n444284_regenie.txt \
    -og - \
    -ofiletype dosage > output/${phenotype}/snps_${chr}_dosage.txt

fi
