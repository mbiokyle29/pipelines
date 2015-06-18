#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for all things tophat

"""
from ruffus import *
from tophat_extras import check_default_args, make_fastq_list
import ruffus.cmdline as cmdline
import subprocess
import logging
import os
import pprint
import re

parser = cmdline.get_argparse(description='Given a directory of NON-paired end reads -- Align them with tophat and generate wigs')

# Program arguments
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--cores", help="Number of cores to run bowtie on", default=10)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file", required=True)
parser.add_argument("--gtf", help="Fullpath to gtf file", required=True)
parser.add_argument("--paired", help="Indicates whether the reads in --dir are paired_end. MUST FOLLOW _1 _2 convention", default=False)

# optional arguments to control turning on and off tasks
parser.add_argument("--wig", help="Whether or not wig files should be generated", type=bool, default=False)
parser.add_argument("--one-codex", help="Whether or not to upload each sample to one codex for metagenomic analysis", default=False)
parser.add_argument("--de", help="Whether or not differential expression should be calculated", type=bool, default=False)
parser.add_argument("--de-conf", help="fullpath to differential expresssion configuration file")


# parse the args
options = parser.parse_args()

# logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO,)
log = logging.getLogger(__name__)
log.info("Starting Bowtie2 Run")

# pre checking
check_default_args(options.cores, options.index, options.output, log)
input_files = make_fastq_list(options.dir, log)

# need this for wig headers
genome = os.path.splitext(os.path.basename(options.index))[0]

@active_if(options.one_codex)
@transform(input_files, suffix(".fastq"), ".fastq")
def upload_to_one_codex(input_file, output_file):

    args = ["onecodex", "upload", input_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("uploading to One Codex failed")
        raise SystemExit

@active_if(options.paired)
@collate(input_files, formatter("([^/]+)_[12].fastq$"), ["{path[0]}/{1[0]}_1.fastq", "{path[0]}/{1[0]}_2.fastq"],)
def collate_files(input_files, output_files):
    log.info("Collating paired fastq files: \n\t{} \n\t{}\n".format(input_files[0], input_files[1]))

# paired alignment
@active_if(options.paired)
@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}-tophat-results/accepted_hits.bam", options)
def tophat_align_paired(input_files, output_file, options):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    args = ["tophat2", "-G", options.gtf,"-p", str(options.cores), "-o", output, options.index, input_files[0], input_files[1]]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("tophat2 failed")
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# unpaired alignment function
@active_if(not options.paired)
@transform(input_files, formatter(), options.output+"{basename[0]}-res/accepted_hits.bam", options)
def tophat_align_unpaired(input_file, output_file, options):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    args = ["tophat2", "-G", options.gtf,"-p", str(options.cores), "-o", output, options.index, input_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("tophat2 failed")
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# both the tophat functions can feed into here
@transform([tophat_align_unpaired, tophat_align_paired], suffix(".bam"),".sorted.bam")
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
@transform(bam_to_bed, suffix(".bed"), ".cov", options.size)
def bed_to_cov(input_file, output_file,  size_file):
    log.info("Converting %s to a genome coverage file", input_file)
    
    command = "genomeCoverageBed -d -i {} -g {} > {}".format( input_file, size_file, output_file)
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

@active_if(options.de)
@merge()

# run the pipeline
cmdline.run(options)