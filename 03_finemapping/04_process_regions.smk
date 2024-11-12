import yaml
import os
import glob
from itertools import product
import json
import sys
import re

#######
# Snakemake workflow to pre-process the regions that are needed for Susie finemapping.
# Per region:
#   list of variants
#   info on variants
#   LD matrix
#######

# Rules that dont need submitting to the cluster
localrules: createMfile

# REGIONS_FILE = "input/testsmall_region.txt" # One small region
REGIONS_FILE = "input/finemapping_regions.txt" # The full thing

info = {}
regions = []
with open(REGIONS_FILE,'r') as infile:
    for line in infile:
        idx, phenotype, region, chrom, start, end = line.strip().split()

        if region in regions:
            continue

        info[region] = {'idx': idx, 'phenotype': phenotype, 'region': region, 'chrom': chrom, 'start': start, 'end': end}
        regions.append(region)

        # Output directory for every region separately
        if not os.path.isdir(f"output/region_data/{region}/"):
            os.mkdir(f"output/region_data/{region}/")

print('Regions:')
print(regions)

print(f"Nr of regions: {len(regions)}")

rule all:
    input:
        expand('output/region_data/{region}/snplist.txt', region=regions),
        expand('output/region_data/{region}/snpdata.z', region=regions),
        expand('output/region_data/{region}/ldmat.ld', region=regions)

rule getSNPlist:
    "Getting list of SNPS in the specific region to extract"
    input:
    output:
        snplist = "output/region_data/{region}/snplist.txt",
    params:
        phenotype = lambda wcs: info[wcs.region]['phenotype'],
        chrom = lambda wcs: info[wcs.region]['chrom'],
        start = lambda wcs: info[wcs.region]['start'],
        end = lambda wcs: info[wcs.region]['end']
    resources:
        mem = "2G",
        runtime = "5m"
    shell:
        """
        zcat output/gwas_collated/gwas_{params.phenotype}_allchr.txt.gz | \
            awk -v chr={params.chrom} -v end={params.end} -v start={params.start} '{{ if (($1==chr && $2 >=start && $2 <= end)) print $3}}' > {output.snplist}
        """

rule getSNPdata:
    "Getting the data on the variants in the region for LDstore"
    input:
    output:
        snpdata = "output/region_data/{region}/snpdata.z"
    params:
        phenotype = lambda wcs: info[wcs.region]['phenotype'],
        chrom = lambda wcs: info[wcs.region]['chrom'],
        start = lambda wcs: info[wcs.region]['start'],
        end = lambda wcs: info[wcs.region]['end']
    resources:
        mem = "5G",
        runtime = "10m"
    shell:
        """
        echo "rsid chromosome position allele1 allele2" > {output.snpdata}

        zcat  output/gwas_collated/gwas_{params.phenotype}_allchr.txt.gz | \
            awk -v chr={params.chrom} -v end={params.end} -v start={params.start} '{{ if (($1==chr && $2 >=start && $2 <= end)) print $0}}' | \
            awk '{{
              if ( $1 < 10) {{
                print $3, "0"$1, $2, $4, $5 ;           # Add leading zero when needed
              }} else if ($1 == 23) {{
                print $3, "X", $2, $4, $5               # Change 23 to X
              }} else {{
                print $3, $1, $2, $4, $5                # Default case
              }}
            }}' >> {output.snpdata}
        """

rule filterBGENvariants:
    "Filtering the original bgen for variants that we need in the region"
    input:
        snplist = "output/region_data/{region}/snplist.txt",
    output:
        fltbgen_variants = "output/region_data/{region}/flt_variants.bgen"
    params:
        chrom = lambda wcs: info[wcs.region]['chrom']
    resources:
        mem = "5G",
        runtime = "15m"
    shell:
        """
        if [ {params.chrom} -eq 23 ]; then
            chr='X'
        else
            chr={params.chrom}
        fi

        bgenix \
            -g /sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c${{chr}}_b0_v3.bgen \
            -incl-rsids {input.snplist} > {output.fltbgen_variants}

        ## Create index for the bgen
        bgenix -g {output.fltbgen_variants} -index -clobber

        """

rule filterBGENsamples:
    "Filter the BGEN (that was already filtered on variants to include) for the samples that we randomly selected"
    input:
        fltbgen_variants = "output/region_data/{region}/flt_variants.bgen"
    output:
        fltbgen_variants_samples = "output/region_data/{region}/flt_variants_samples.bgen"
    params:
        samplefile = "output/randomsamples.incl",
        chrom = lambda wcs: info[wcs.region]['chrom']
    resources:
        mem = "15G",
        runtime = "5h"
    shell:
        """
        if [ {params.chrom} -eq 23 ]; then
            chr='X'
        else
            chr={params.chrom}
        fi

        /programs/qctool/build/release/qctool_v2.0.7 \
            -g {input.fltbgen_variants} \
            -s /sc-resources/ukb/data/bulk/genetic/imputed/bgen_files_44448/ukb22828_c${{chr}}_b0_v3.sample \
            -incl-samples {params.samplefile} -og {output.fltbgen_variants_samples}

        programs/bgen/build/apps/bgenix -g {output.fltbgen_variants_samples} -index -clobber
        """


rule createMfile:
    "Create the masterfile for LDStore"
    input:
        snpdata =  "output/region_data/{region}/snpdata.z",
        bgen = "output/region_data/{region}/flt_variants_samples.bgen"
    output:
        mfile = "output/region_data/{region}/mfile.z"
    params:
        chrom = lambda wcs: info[wcs.region]['chrom'],
        ldmat = "output/region_data/{region}/ldmat.ld",
        samplefile = "output/randomsamples.incl",
        bgi = "output/region_data/{region}/flt_variants_samples.bgen.bgi"
    shell:
        """
        if [ {params.chrom} -eq 23 ]; then
            n=49948
        else
            n=50000
        fi

        echo "z;bgen;bgi;ld;incl;n_samples" > {output.mfile}
        echo "{input.snpdata};{input.bgen};{params.bgi};{params.ldmat};{params.samplefile};$n" >> {output.mfile}

        """

rule calcLD:
    "Finally perform the actual LD calculations using LDStore"
    input:
        mfile = "output/region_data/{region}/mfile.z"
    output:
        ldmat = "output/region_data/{region}/ldmat.ld"
    params:
    resources:
        cpus_per_task = 15,
        mem = '75G',
        runtime = '12h'
    shell:
        """
        programs/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
            --in-files {input.mfile} \
            --write-text \
            --n-threads 15 \
            --read-only-bgen
            """
