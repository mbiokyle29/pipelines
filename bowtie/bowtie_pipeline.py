#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for simple bowtie alignment

"""
from ruffus import *
from bowtie_extras import BowtieExtras
import ruffus.cmdline as cmdline
import subprocess
import logging
import os
import pprint
import re
import time

parser = cmdline.get_argparse(description='Given a directory of NON-paired end reads -- Align them with bowtie')

# Program arguments  -- Most go straight to bowtie
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True, nargs='+')
parser.add_argument("--cores", help="Number of cores to run bowtie on", default=10)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file")

# optional arguments to control turning on and off tasks
parser.add_argument("--wig", help="Whether or not wig files should be generated", type=bool, default=False)

# parse the args
options = parser.parse_args()

# Kenny loggins
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

# file log
time_stamp = str(time.time()).replace(".","")
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("bowtie_pipeline",time_stamp,"log"))
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(log_formatter)

# console log
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(log_formatter)

# set it all up
log.addHandler(file_handler)
log.addHandler(stream_handler)
log.info("Starting Bowtie Run")

extras = BowtieExtras(log)

# alignment stats file (kind of not pipelined)
for parameter in ["dir","cores","index","output"]:
    log.info("{}: {}".format(parameter,getattr(options,parameter)))

# pre checking
extras.check_default_args(options.cores, options.index, options.output)
input_files = extras.make_fastq_list(options.dir)

# storing MRM array
stats_file = os.path.join(options.output, "bowtieReadCounts-{}.stats".format(time_stamp))

# need this for wig headers
genome = os.path.splitext(os.path.basename(options.index))[0]

@transform(input_files, suffix(".fastq"), ".sam", options, stats_file, extras)
def align_with_bowtie(input_file, output_file, options, stats_file, extras):
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
    extras.record_bowtie_output(output, input_file, stats_file)

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
@transform(bam_to_bed, suffix(".bed"), ".bg", options.size, extras)
def bed_to_bg(input_file, output_file,  size_file, extras):
    log.info("Converting %s to a genome coverage file", input_file)
    
    base = os.path.splitext(os.path.basename(input_file))[0]
    mill_reads = extras.get_mill_reads(base)
    scale = 1 / mill_reads

    command = "genomeCoverageBed -scale {} -bg -split -i {} -g {} > {}".format(scale, input_file, size_file, output_file)
    if subprocess.call(command, shell=True):
        log.warn("bed to coverage conversion of %s failed, exiting", input_file)
        extras.report_error("bed_to_bg","bed to bg conversion of {} failed".format(input_file))
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@active_if(options.wig)
@transform(bed_to_bg, formatter(), options.output+"{basename[0]}.bw", options.size, options.output)
def bg_to_bw(input_file, output_file, size_file, output):
    log.info("Creating bigwig file from bg:  %s", input_file)
    command = "bedGraphToBigWig {} {} {}".format( input_file, size_file, output_file)
    
    if subprocess.call(command, shell=True):
        log.warn("bg to bw conversion of %s failed, exiting", input_file)
        extras.report_error("bg_to_bw","bg to bw conversion of {} failed".format(input_file))
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)


# call out to external bwtools here
@active_if(options.wig)
@merge(bg_to_bw, os.path.join(options.output,"bigWigStats-"+time_stamp+".out"))
def bw_stats(input_files, output_file):
        
    # we are going to call bwtool summary and bwtool distribution
    # have to explicitly send stdout stuff like that
    # what a program
    summary = "bwtool summary 10000 -header -with-sum {} /dev/stdout"
    dist    = "bwtool distribution {} /dev/stdout"
    

    for input_file in input_files:
        log.info("Running bigwig stats on {}".format(input_file))
        with open(output_file, "a+") as stats:
            for command in [summary, dist]:
                
                    try:
                        output = subprocess.check_output(command.format(os.path.abspath(input_file)).split())
                        if command.startswith("bwtool summary"):
                            stats.write("#### bwtool summary for {}\n".format(input_file))
                            stats.write(output)
                            stats.write("####\n")
                        
                        # filter zeros out
                        else:
                            output = output.rstrip()
                            output_clean = [line for line in output.split("\n") if line.split('\t')[1] != '0']
                            stats.write("#### bwtool distribution for {}\n".format(input_file))
                            stats.write("depth\tcount\n")
                            stats.write("\n".join(output_clean))
                            stats.write("\n####\n")

                    except subprocess.CalledProcessError:
                        log.warn("{} failed running on {}".format(command, input_file))
                        raise SystemExit
            stats.write("\n\n")
# run the pipelined
cmdline.run(options)