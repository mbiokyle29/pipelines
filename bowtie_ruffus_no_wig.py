#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for simple bowtie alignment

"""
from ruffus import *
from bowtie_extras import check_default_args, make_fastq_list, record_bowtie_output
import ruffus.cmdline as cmdline
import subprocess
import logging
import os
import pprint
import re

parser = cmdline.get_argparse(description='Given a directory of NON-paired end reads -- Align them with bowtie')

# Program arguments  -- Most go straight to bowtie
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--cores", help="Number of cores to run bowtie on", default=10)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file", required=True)

# parse the args
options = parser.parse_args()

# logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO,)
log = logging.getLogger(__name__)
log.info("Starting Bowtie2 Run")

# alignment stats file (kind of not pipelined)
for parameter in ["dir","cores","index","output"]:
    log.info("{}: {}".format(parameter,getattr(options,parameter)))

# pre checking
check_default_args(options.cores, options.index, options.output, log)
input_files = make_fastq_list(options.dir, log)

# storing MRM array
mrm = {}
stats_file = os.path.join(options.output, "bowtie-pipeline-no-wig.stats")

# need this for wig headers
genome = os.path.splitext(os.path.basename(options.index))[0]


@transform(input_files, suffix(".fastq"),".sam", options, stats_file)
def align_with_bowtie(input_file, output_file, options, stats):
    log.info("Running bowtie2 on %s", input_file)
    
    # use poen explicitly to cpature STDERR, check still
    args = ["bowtie2", "-t", "--no-unal", "-p", str(options.cores), "-x", options.index, input_file, "-S", output_file]
    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("Bowtie failed")
        raise SystemExit

    # print output
    log.info("Bowtie output:")
    log.info(output)

    # pass along to be saved
    mrms = record_bowtie_output(output, input_file, stats_file, log)
    base = os.path.splitext(os.path.basename(input_file))[0]
    mrm[base] = mrms

#
@transform(align_with_bowtie, suffix(".sam"),".bam")
def sam_to_bam(input_file, output_file):

    log.info("Converting %s to bam", input_file)
    if subprocess.call(["samtools", "view", "-S", "-b", "-o", output_file, input_file]):
        log.warn("sam to bam convertion of %s failed, exiting", input_file)
        raise SystemExit

    # clean up
    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@transform(sam_to_bam, suffix(".bam"),".sorted.bam")
def sort_bam(input_file, output_file):
    log.info("Sorting %s ", input_file)

    # hacky
    output_file = re.sub(r"\.bam", "", output_file)

    if subprocess.call(["samtools-rs", "rocksort", "-@", "8", "-m", "16G", input_file, output_file]):
        log.warn("bam sorting %s failed, exiting", input_file)
        raise SystemExit
    
    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

# run the pipeline
cmdline.run(options)