import yaml
import os
import glob
from itertools import product
import json
import sys
import re

localrules: createMETALscript

# Snakemake workflow to ru meta analysis across the three ancestries

# Contains one single column of the exact phenotype names that were mapped.
# Can be created by using the phenotype files used as input for REGENIE step 2
# e.g.: zcat input/phenotypes.txt.gz | head -n1 | tr '\t' '\n' | grep 'FID\|IID' -v > input/phenotypes_mapped.txt

PHENOTYPE_FILE = "input/phenotypes.txt"
# CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

PHENOTYPES = []

# if not os.path.isdir("output/meta_analyses/data/"):
#     os.mkdir("output/meta_analyses/data/")
#     os.mkdir("output/meta_analyses/data/eur/")
#     os.mkdir("output/meta_analyses/data/csa/")
#     os.mkdir("output/meta_analyses/data/afr")

if not os.path.isdir("output/meta_analyses/scripts/"):
    os.mkdir("output/meta_analyses/scripts/")

if not os.path.isdir("output/meta_analyses/meta/"):
    os.mkdir("output/meta_analyses/meta/")

if not os.path.isdir("output/meta_analyses/clump/"):
    os.mkdir("output/meta_analyses/clump/")

    if not os.path.isdir("output/meta_analyses/collated_clumps/"):
        os.mkdir("output/meta_analyses/collated_clumps/")

with open(PHENOTYPE_FILE,'r') as infile:
    for line in infile:
        phenotype = line.strip()

        PHENOTYPES.append(phenotype)

rule all:
    input:
        expand("output/meta_analyses/meta_flt/meta_flt_{phenotype}.txt.gz", phenotype=PHENOTYPES),
        expand("output/regions_interest/regions_{phenotype}.txt", phenotype=PHENOTYPES),
        expand("output/sentinels/sentinels_{phenotype}.txt", phenotype=PHENOTYPES)

rule createMETALscript:
    input:
        eur = 'output/gwas_collated/gwas_{phenotype}_allchr.txt.gz',
        csa = '01_gwas_pop_csa/output/gwas_collated/gwas_{phenotype}_allchr.txt.gz',
        afr = '02_gwas_pop_afr/output/gwas_collated/gwas_{phenotype}_allchr.txt.gz'
    output: outfile = "output/meta_analyses/scripts/script_{phenotype}.txt"
    params: prefix = "output/meta_analyses/meta/meta_{phenotype}"
    shell:
        """
        # Create the METAL scripts to be used

        # Inverse variance-weighted method that uses effect sizes and standard errors


        echo "
        SCHEME   STDERR
        STDERRLABEL   SE
        TRACKPOSITIONS ON
        LOGPVALUE ON

        MARKER   ID
        WEIGHT   N
        ALLELE   ALLELE0 ALLELE1
        FREQ     A1FREQ
        EFFECT   BETA
        STDERR   SE
        CHROMOSOME  CHROM
        POSITION GENPOS

        PROCESS {input.eur}
        PROCESS {input.csa}
        PROCESS {input.afr}

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

rule filterMETALouts:
    input:
        metalin = "output/meta_analyses/meta/meta_{phenotype}1.txt",
        snpfile = 'output/common_variants.txt'
    params:
        tmpfile = 'output/meta_analyses/meta_flt/meta_flt_{phenotype}.txt'
    output:
        metalout = 'output/meta_analyses/meta_flt/meta_flt_{phenotype}.txt.gz'
    resources:
        mem='2G',
        runtime='25m'
    shell:
        """
        grep -Fw -f {input.snpfile} {input.metalin} > {params.tmpfile}
        gzip {params.tmpfile}
        """

rule get_regions_interest:
    input: metalin = 'output/meta_analyses/meta_flt/meta_flt_{phenotype}.txt.gz'
    output: regions_of_interest = "output/regions_interest/regions_{phenotype}.txt"
    params:
        gw_threshold = -9.7
    resources:
        mem='2G',
        runtime='30m'
    shell:
        """
        ## get all signals of interest (treat extended MHC region as one)
        # Adapted for METAL input

        # THe logged P-value that METAL reports is the log10P and not the -log10P
        # Therefore, filter for results that are < -9.7, equal to > 9.7
        # Also no need to skip the first line, no header

        zcat {input.metalin} | \
        awk '{{if (($8+0) < {params.gw_threshold}) print $0}}'| \
        awk '{{print $1,$2,$2}}' OFS="\t"  | \
        sort -k1,1 -k2,2n | \
        awk '{{if($1+0.0==6+0 && $2+0 >= 25500000+500000 && $2+0 <= 34000000-500000) {{print $1,25500000,34000000}} \
                 else if($1+0.0==6+0 && $2+0 >= 25500000 && $2+0 <= 25500000+500000) {{print $1,$2-500000,34000000}} \
                 else if($1+0.0==6+0 && $2+0 >= 34000000-500000 && $2+0 <= 34000000) {{print $1,25500000,$2+500000}} \
                 else if($2-500000 >= 0) {{print $1,$2-500000,$2+500000}} \
                 else {{print $1,0,$2+500000}} }}' OFS="\t" | \
        bedtools merge -i stdin > {output.regions_of_interest}
        """



rule define_regional_sentinels:
    input:
        metalin = 'output/meta_analyses/meta_flt/meta_flt_{phenotype}.txt.gz',
        regions = "output/regions_interest/regions_{phenotype}.txt"
    output:
        sentinels = "output/sentinels/sentinels_{phenotype}.txt"
    params:
        gw_threshold = -9.7
    resources:
        mem='20G',
        runtime='1h'
    shell:
        """
            # Get the strongest signal for each of the regions
            # do in R for easy way

            # Variables needed for singularity
            R_CONTAINER='all_inclusive_rstudio_4.3.2.sif'
            R_SCRIPT='scripts/01_getsentinels.R'
            BIND_DIR=""

            singularity exec \
              --bind $BIND_DIR \
              $R_CONTAINER Rscript $R_SCRIPT {input.metalin} {input.regions} {output.sentinels}
        """
