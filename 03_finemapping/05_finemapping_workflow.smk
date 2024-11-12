import yaml
import os
import glob
from itertools import product
import json
import sys
import re

# Snakemake workflow to perform fine mapping with SusieR
# Requires regional data to be ready (e.g. 04_process_regions.smk)

# REGIONS_FILE = "input/testsmall.txt" # Test cases
REGIONS_FILE = "input/finemapping_regions.txt" # The full thing

if not os.path.isdir(f"output/finemapping/"):
    os.mkdir(f"output/finemapping/")

info = {}
indices = []
with open(REGIONS_FILE,'r') as infile:
    for line in infile:
        idx, phenotype, region, chrom, start, end = line.strip().split()

        info[idx] = {'idx': idx, 'phenotype': phenotype, 'region': region, 'chrom': chrom, 'start': start, 'end': end}
        indices.append(idx)

        # Output directory for every region separately
        if not os.path.isdir(f"output/finemapping/{idx}/"):
            os.mkdir(f"output/finemapping/{idx}/")

rule all:
    input:
        expand('output/finemapping/{idx}/finemapped_results.txt.gz', idx = indices)

def get_mem(wildcards, attempt):
    # print(f"This is attempt: {attempt}")
    if attempt == 1:
        return '75G'
    elif attempt == 2:
        return '150G'
    elif attempt == 3:
        return '250G'

rule runFinemapping:
    input:
    output: "output/finemapping/{idx}/finemapped_results.txt.gz"
    params:
        index = lambda wcs: info[wcs.idx]['idx'],
        phenotype = lambda wcs: info[wcs.idx]['phenotype'],
        region = lambda wcs: info[wcs.idx]['region'],
        chrom = lambda wcs: info[wcs.idx]['chrom'],
        start = lambda wcs: info[wcs.idx]['start'],
        end = lambda wcs: info[wcs.idx]['end']
    resources:
        mem = get_mem,
        time = "36:00:00",
        cpus_per_task = 5
    shell:
        """
            # Variables needed for singularity
            R_CONTAINER='programs/all-inclusive-rstudio-apptainer/sif/all_inclusive_rstudio_4.3.2.sif'
            R_SCRIPT='scripts/05_run_susie.R'
            BIND_DIR=""

            singularity exec \
              --bind $BIND_DIR \
              $R_CONTAINER Rscript $R_SCRIPT {params.index} {params.phenotype} {params.region} {params.chrom} {params.start} {params.end}
        """
