#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for simple bowtie alignment

"""
from ruffus import *
from bowtie_extras import check_default_args, record_bowtie_output, make_fastq_list
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

# optional arguments to control turning on and off tasks
parser.add_argument("--wig", help="Whether or not wig files should be generated", type=bool, default=False)

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
stats_file = os.path.join(options.output, "bowtie-pipeline.stats")

# need this for wig headers
genome = os.path.splitext(os.path.basename(options.index))[0]

@transform(input_files, suffix(".fastq"), ".sam", options, stats_file)
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

@active_if(options.wig)
@transform(sort_bam, suffix(".sorted.bam"), ".bed", options.output)
def bam_to_bed(input_file, output_file, output):
    log.info("Converting %s to a bed file", input_file)
    if subprocess.call("bamToBed -i {} > {}".format(input_file, output_file), shell=True):
        log.warn("bam to bed conversion of %s failed, exiting", input_file)
        raise SystemExit

    # now we can move sorted bam to output
    file_name = os.path.basename(input_file)
    new_name = os.path.join(output, file_name)
    os.rename(input_file, new_name)

@active_if(options.wig)
@transform(bam_to_bed, suffix(".bed"), ".cov", mrm, options.size)
def bed_to_cov(input_file, output_file, mrms, size_file):
    log.info("Converting %s to a genome coverage file", input_file)
    
    # check
    base = os.path.splitext(os.path.basename(input_file))[0]
    scale = mrms[base]

    command = "genomeCoverageBed -scale {} -d -i {} -g {} > {}".format(scale, input_file, size_file, output_file)
    if subprocess.call(command, shell=True):
        log.warn("bed to coverage conversion of %s failed, exiting", input_file)
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@active_if(options.wig)
@transform(bed_to_cov, suffix(".cov"), ".wig", genome, options.output)
def cov_to_wig(input_file, output_file, genome, output):
    log.info("Creating wig file from coverage bed %s", input_file)

    output_stream = open(output_file,"w+")

    # write the header
    base = os.path.splitext(os.path.basename(input_file))[0]
    desc = "{} aligned to {} with bowtie ruffus pipeline".format(base, genome)
    header = "track type=wiggle_0 name=\"{}\" description=\"{}\"\n".format(base,desc)
    track = "fixedStep chrom={} start=1 step=1\n".format(genome)

    output_stream.write(header)
    output_stream.write(track)

    with open(input_file,"r") as input:
        for line in input:
            output_stream.write(line.split("\t")[2])

    input.close()
    output_stream.close()

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

    # move the wig
    file_name = os.path.basename(output_file)
    new_name = os.path.join(output, file_name)
    os.rename(output_file, new_name)

# run the pipeline
cmdline.run(options)