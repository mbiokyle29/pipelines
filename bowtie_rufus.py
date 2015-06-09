#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for simple bowtie alignment

"""
from ruffus import *
from bowtie_extras import check_default_args, make_fastq_list
import ruffus.cmdline as cmdline
import subprocess
import logging
import os
import pprint

parser = cmdline.get_argparse(description='Given a directory of NON-paired end reads -- Align them with bowtie')

# Program arguments  -- Most go straight to bowtie
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--cores", help="Number of cores to run bowtie on", default=10)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")

# parse the args
options = parser.parse_args()

# logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
log = logging.getLogger(__name__)
log.info("Starting Bowtie2 Run")

for parameter in ["dir","cores","index","output"]:
    log.info("{}: {}".format(parameter,getattr(options,parameter)))

# pre checking
check_default_args(options.cores, options.index, options.output, log)
input_files = make_fastq_list(options.dir, log)

# Step 1 uncompress (bowtie cant take gz)
@transform(input_files, suffix(".fastq.gz"),".fastq")
def gunzip(input_file, output_file):
    
    log.info("gunzipping %s", input_file)   

    # 0 == all good, and bool false
    if subprocess.call(["gunzip",input_file]):
        log.warn("gunzipping %s failed, exiting", input_file)
        raise SystemExit

# Step 2 align
@transform(gunzip, suffix(".fastq"),".sam", options)
def align_with_bowtie(input_file, output_file, options):
    log.info("Running bowtie2 on %s", input_file)
    if subprocess.call(["bowtie2", "-t", "--no-unal", "-p", str(options.cores), "-x", options.index, input_file, "-S", output_file]):
        log.warn("bowtie alignment of %s failed, exiting", input_file)
        raise SystemExit

# step 3 BAM
@transform(align_with_bowtie, suffix(".sam"),".bam")
def sam_to_bam(input_file, output_file):
    log.info("Converting %s to bam", input_file)
    if subprocess.call(["samtools", "view", "-S", "-b", "-o", output_file, input_file]):
        log.warn("sam to bam convertion of %s failed, exiting", input_file)
        raise SystemExit

# step 4
@transform(sam_to_bam, suffix(".bam"),".sorted")
def sort_bam(input_file, output_file):
    log.info("Sorting %s ", input_file)
    if subprocess.call(["samtools-rs", "rocksort", "-@", "8", "-m", "16G", input_file, output_file]):
        log.warn("bam sorting %s failed, exiting", input_file)
        raise SystemExit

# run the pipeline
cmdline.run(options)