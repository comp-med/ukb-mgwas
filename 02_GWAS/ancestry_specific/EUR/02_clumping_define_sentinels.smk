import yaml
import os
import glob
from itertools import product
import json
import sys
import re

# Snakemake workflow to perform:
#   1. Regional clumping
#   2. Fine-mapping

# Contains one single column of the exact phenotype names that were mapped.
# Can be created by using the phenotype files used as input for REGENIE step 2
# e.g.: zcat input/phenotypes.txt.gz | head -n1 | tr '\t' '\n' | grep 'FID\|IID' -v > input/phenotypes_mapped.txt

GWAS_LOCATION="output/gwas_collated/"
COMBINATIONS_FILE = "input/phenotypes_mapped.txt"
GWSIG_THRESHOLD = 9.7 # -log10(5e-8 / 249)

PHENOTYPES = []

with open(COMBINATIONS_FILE,'r') as infile:
    for line in infile:
        phenotype = line.strip()
        PHENOTYPES.append(phenotype)

rule all:
    input:
        regions = expand("output/regions_interest/regions_{phenotype}.txt", phenotype=PHENOTYPES),
        sentinels = expand("output/sentinels/sentinels_{phenotype}.txt", phenotype=PHENOTYPES)


rule get_regions_interest:
    input: summarystats = GWAS_LOCATION + "gwas_{phenotype}_allchr.txt.gz"
    output: regions_of_interest = "output/regions_interest/regions_{phenotype}.txt"
    params:
        gw_threshold = GWSIG_THRESHOLD
    resources:
        mem='10G',
        runtime='10m'
    shell:
        """
        ## get all signals of interest (treat extended MHC region as one)

        zcat {input.summarystats} | \
        awk '{{if(NR == 1) print $0; else if (NR > 1 && ($13+0) > {params.gw_threshold}) print $0}}'| \
        sed '1d'  | \
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
        summarystats = GWAS_LOCATION + "gwas_{phenotype}_allchr.txt.gz",
        regions = "output/regions_interest/regions_{phenotype}.txt"
    output:
        sentinels = "output/sentinels/sentinels_{phenotype}.txt"
    params:
        gw_threshold = GWSIG_THRESHOLD
    resources:
        mem='10G',
        runtime='3h'
    shell:
        """
            # Get the strongest signal for each of the regions

            while read chr start end
            do

            zcat {input.summarystats} | \
            awk '{{if(NR==1) print $0; else if(NR>1 && ($13+0.0) > {params.gw_threshold}) print $0}}' OFS=" " | \
            sed '1d' | \
            awk -v chr=${{chr}} -v start=${{start}} -v end=${{end}} '{{if(($1+0)==(chr+0) && ($2+0) >= (start+0) && ($2+0) <= (end+0)) print}}' | \
            awk -v start=${{start}} -v end=${{end}} -v max=0 '{{if(sqrt(($10/$11)*($10/$11))>=max){{want=$0; max=sqrt(($10/$11)*($10/$11))}}}} END{{print want,start,end}}' OFS=" " ;

            done < {input.regions} > {output.sentinels}
        """
