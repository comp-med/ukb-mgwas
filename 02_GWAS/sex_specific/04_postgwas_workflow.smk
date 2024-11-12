import yaml
import os
import glob
from itertools import product
import json
import sys
import re

# Snakemake workflow to collate the output from REGENIE, filter variants of low quality (INFO and MAF) and run LDSC on each collated summary statistic

# Contains one single column of the exact phenotype names that were mapped.
# Can be created by using the phenotype files used as input for REGENIE step 2
# e.g.: zcat input/phenotypes.txt.gz | head -n1 | tr '\t' '\n' | grep 'FID\|IID' -v > input/phenotypes_mapped.txt

GWSIG_THRESHOLD = 9.7 # -log10(5e-8 / 249)

COMBINATIONS_FILE = "input/phenotypes_mapped.txt"
CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']
SEXES=['males', 'females']

PHENOTYPES = []

with open(COMBINATIONS_FILE,'r') as infile:
    for line in infile:
        phenotype = line.strip()
        PHENOTYPES.append(phenotype)

rule all:
    input:
        expand("output/gwas_collated/{sex}/gwas_{phenotype}_allchr.txt.gz", phenotype=PHENOTYPES, sex=SEXES),
        # expand("output/ld_score_regression/{phenotype}.log", phenotype=PHENOTYPES),
        # expand("output/ld_score_regression/{phenotype}_h2.log", phenotype=PHENOTYPES),
        expand("output/regions_interest/{sex}/regions_{phenotype}.txt", phenotype=PHENOTYPES, sex=SEXES),
        expand("output/sentinels/{sex}/sentinels_{phenotype}.txt", phenotype=PHENOTYPES, sex=SEXES)

rule collate:
    input: infiles = expand("output/gwas/{{sex}}/chr{chr}_{{phenotype}}.regenie.gz", chr=CHROMOSOMES)
    output: outfile = "output/gwas_collated/{sex}/gwas_{phenotype}_allchr.txt.gz"
    params: tmpfile = "output/gwas_collated/{sex}/gwas_{phenotype}_allchr.txt"
    resources:
        mem='3G',
        runtime='15m'
    shell:
        """
        # Create header for the output file
        echo -e 'CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA' > {params.tmpfile}

        # Remove the header per file and filter variants
        # MAF: >0.5%
        # INFO: >0.5

        zcat {input.infiles} | \
        grep -v 'CHROM' | \
        awk '{{ if( ($7 >= 0.5) &&
                    ($6 >= 0.005 && $6 <= 0.995)) print $0 }}' >> {params.tmpfile}

        gzip {params.tmpfile}
        """

rule ldsc:
    input: summarystats = "output/gwas_collated/gwas_{phenotype}_allchr.txt.gz"
    output:
        outfile = "output/ld_score_regression/{phenotype}.log",
        outfile_h2 = "output/ld_score_regression/{phenotype}_h2.log"
    params: tmpfile = "tmp/{phenotype}.tmp"
    conda: "ldsc"
    resources:
        mem='25G',
        runtime='15m'
    shell:
        """

        # Reformat the REGENIE output to something temporary
        zcat {input.summarystats} | \
        awk '{{if(NR==1) print $3"\t"$1"\t"$2"\t"$5"\t"$4"\t"$6"\t"$10"\t"$11"\t"$13"\t"$8; else print $3"\t"$1"\t"$2"\t"$5"\t"$4"\t"$6"\t"$10"\t"$11"\t"10^-$13"\t"$8}}' - > {params.tmpfile}

        # Run LD-score regression
        ldsc/munge_sumstats.py \
            --sumstats {params.tmpfile} \
            --snp ID \
            --a1 ALLELE1 \
            --a2 ALLELE0 \
            --p LOG10P \
            --out "output/ld_score_regression/{wildcards.phenotype}" \
            --chunksize 500000 \
            --merge-alleles programs/ldsc/w_hm3.snplist

        programs/ldsc/ldsc.py \
            --h2 output/ld_score_regression/{wildcards.phenotype}.sumstats.gz \
            --ref-ld-chr programs/ldsc/eur_w_ld_chr/ \
            --w-ld-chr programs/ldsc/eur_w_ld_chr/ \
            --out "output/ld_score_regression/{wildcards.phenotype}_h2"


        # Get rid of tmp files
        rm {params.tmpfile}
        rm output/ld_score_regression/{wildcards.phenotype}.sumstats.gz
        """

rule get_regions_interest:
    input: summarystats = "output/gwas_collated/{sex}/gwas_{phenotype}_allchr.txt.gz"
    output: regions_of_interest = "output/regions_interest/{sex}/regions_{phenotype}.txt"
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
        summarystats = "output/gwas_collated/{sex}/gwas_{phenotype}_allchr.txt.gz",
        regions = "output/regions_interest/{sex}/regions_{phenotype}.txt"
    output:
        sentinels = "output/sentinels/{sex}/sentinels_{phenotype}.txt"
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
