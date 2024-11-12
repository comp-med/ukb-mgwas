import yaml
import os
import glob
from itertools import product
import json
import sys
import re

localrules: createMETALscript, collate

# Snakemake workflow to collate the output from REGENIE, filter variants of low quality (INFO and MAF) and run LDSC on each collated summary statistic

# Contains one single column of the exact phenotype names that were mapped.
# Can be created by using the phenotype files used as input for REGENIE step 2
# e.g.: zcat input/phenotypes.txt.gz | head -n1 | tr '\t' '\n' | grep 'FID\|IID' -v > input/phenotypes_mapped.txt

COMBINATIONS_FILE = "input/phenotypes_mapped.txt"
CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']
SEXES=['males', 'females']

PHENOTYPES = []

with open(COMBINATIONS_FILE,'r') as infile:
    for line in infile:
        phenotype = line.strip()

        # Output directory for every region separately
        if not os.path.isdir(f"output/meta_analyses/clump/{phenotype}/"):
            os.mkdir(f"output/meta_analyses/clump/{phenotype}/")


        PHENOTYPES.append(phenotype)

# Test run
# PHENOTYPES = ['Ala', 'Glucose', 'M_VLDL_P']

rule all:
    input:
        expand("output/meta_analyses/meta/meta_{phenotype}1.txt", phenotype=PHENOTYPES),
        expand("output/meta_analyses/clump/{phenotype}/clump_chr{chr}.txt", phenotype=PHENOTYPES, chr=CHROMOSOMES),
        expand("output/meta_analyses/collated_clumps/{phenotype}_clumps.txt", phenotype=PHENOTYPES)


rule convert:
    input:
    output:
        males = "output/meta_analyses/data/males/{phenotype}.txt",
        females = "output/meta_analyses/data/females/{phenotype}.txt"
    resources:
        mem='25G',
        runtime='30m'
    shell:
        """
        # Variables needed for singularity
        R_CONTAINER='all-inclusive-rstudio-apptainer/sif/all_inclusive_rstudio_4.3.2.sif'
        R_SCRIPT='scripts/05_process_summarystats.R'
        BIND_DIR=""

        singularity exec \
          --bind $BIND_DIR \
          $R_CONTAINER Rscript $R_SCRIPT {wildcards.phenotype}
        """

rule createMETALscript:
    input:
        males = "output/meta_analyses/data/males/{phenotype}.txt",
        females = "output/meta_analyses/data/females/{phenotype}.txt"
    output: outfile = "output/meta_analyses/scripts/script_{phenotype}.txt"
    params: prefix = "output/meta_analyses/meta/meta_{phenotype}"
    shell:
        """
        # Create the METAL scripts to be used

        # Inverse variance-weighted method that uses effect sizes and standard errors


        echo "
        SCHEME   STDERR
        STDERRLABEL   SE


        MARKER   ID
        WEIGHT   N
        ALLELE   ALLELE0 ALLELE1
        FREQ     A1FREQ
        EFFECT   BETA
        STDERR   SE
        PVAL     P

        PROCESS {input.males}

        MARKER   ID
        WEIGHT   N
        ALLELE   ALLELE0 ALLELE1
        FREQ     A1FREQ
        EFFECT   BETA
        STDERR   SE
        PVAL     P

        PROCESS {input.females}

        OUTFILE {params.prefix}  .txt
        ANALYZE  HETEROGENEITY

        QUIT
        " > {output.outfile}
        """

rule METAL:
    input: metalscript = "output/meta_analyses/scripts/script_{phenotype}.txt"
    output: metalout = "output/meta_analyses/meta/meta_{phenotype}1.txt"
    resources:
        mem='30G',
        runtime='1h'
    shell:
        """
        metal < {input.metalscript}
        """

def get_mem(wildcards, attempt):
    print(f"This is attempt: {attempt}")
    if attempt == 1:
        return '25G'
    elif attempt == 2:
        return '35G'
    elif attempt == 3:
        return '50G'

rule clump_singlechr:
    input: metalout = "output/meta_analyses/meta/meta_{phenotype}1.txt"
    output: out = "output/meta_analyses/clump/{phenotype}/clump_chr{chr}.txt"
    params: varlist = "output/meta_analyses/variants/varlist_{phenotype}.txt"
    resources:
        mem=get_mem,
        runtime='1h'
    shell:
        """
            plink2=programs/plink/plink2_linux_amd_avx2_20240105/plink2

            # Extract list of unique variants that we consider

            # Clumping here on only P-value and not LD.
            # Noticed that some peaks of the heterogeneity Pvalues can be quite wide, this will always lead to spurious SNPs that are
            # for some reason not in LD with the index SNP and therefore form a clump of themselves.
            # For the paper, this is not what we're after and therefore set a (very) broad boundary.
            # THis will pick up only the very largest of signals and smush together any secondary signal whatsoever.
            $plink2 \
                --pfile plink2_files/files/plink2_chr{wildcards.chr} \
                --rm-dup force-first \
                --memory 20000 \
                --clump {input.metalout} \
                --clump-field HetPVal \
                --clump-snp-field MarkerName \
                --clump-p1 5e-8 \
                --clump-kb 2500 \
                --clump-r2 0 \
                --clump-p2 0.0001 \
                --out output/meta_analyses/clump/clump_{wildcards.phenotype}_chr{wildcards.chr}


            # If there is no sig. clump, Plink does not produce an output file.
            # Therefore, touch an empty file to satisfy Snakemake, and stuff these away in separate directories

            touch output/meta_analyses/clump/{wildcards.phenotype}/clump_chr{wildcards.chr}.txt

        """

rule collate:
    input: out = expand("output/meta_analyses/clump/{{phenotype}}/clump_chr{chr}.txt", chr=CHROMOSOMES)
    # out = "output/meta_analyses/clump/{phenotype}/clump_chr{chr}.txt"
    output: file = "output/meta_analyses/collated_clumps/{phenotype}_clumps.txt"
    shell:
        """
            # Get the header
            echo -e "CHROM\tPOS\tID\tP\tTOTAL\tNONSIG\tS0.05\tS0.01\tS0.001\tS0.0001\tSP2" > {output.file}

            cat output/meta_analyses/clump/clump_{wildcards.phenotype}_chr*.clumps |
                grep -v 'CHROM' >> {output.file}
        """
