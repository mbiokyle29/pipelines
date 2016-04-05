#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for all things tophat

"""

# ruffus imports
from ruffus import *
import ruffus.cmdline as cmdline

# custom functions
from tophat_extras import TophatExtras

# system imports
import subprocess, logging, os, re, time

# :) so i never have to touch excel
import pandas as pd

# EMAIL
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

parser = cmdline.get_argparse(description='This pipeline provides a number of funtionalities for working with RNAseq data')

# Program arguments
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--cores", help="Number of cores to run bowtie on", default='10')
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file", default="/data/refs/hg19/hg19.sizes")
parser.add_argument("--gtf", help="Fullpath to gtf file", default="/data/refs/hg19/hg19.gtf")
parser.add_argument("--paired", help="Indicates whether the reads in --dir are paired_end. MUST FOLLOW _1 _2 convention", default=False)

# optional arguments to control turning on and off tasks
# trimming
parser.add_argument("--trim", help="Whether or not to trim the fastq reads", type=bool, default=False)
parser.add_argument("--trim-val", help="Value to trim by (from the end of read)", default=50)

# reporting / meta analysis
parser.add_argument("--bw", help="Whether or not big wig files should be generated", type=bool, default=False)
parser.add_argument("--one-codex", help="Whether or not to upload each sample to one codex for metagenomic analysis", default=False)

# reporting
parser.add_argument("--emails", help="Emails to send DE results too", default="mbio.kyle@gmail.com", nargs="+")

# parse the args
options = parser.parse_args()

# package the emails into an array if just one
if options.emails == "mbio.kyle@gmail.com":
    options.emails = [options.emails]

# Kenny loggins
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

# file log
time_stamp = str(time.time()).replace(".","")
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("tophat_pipeline",time_stamp,"log"))
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
log.info("Starting Tophat")

# pre checking
# create extas instance
extras = TophatExtras(log)

extras.check_default_args(options.cores, options.index, options.output)
input_files = extras.make_fastq_list(options.dir)
genome = os.path.splitext(os.path.basename(options.index))[0]

@active_if(options.paired)
@collate(input_files, formatter("([^/]+)_[12].fastq$"), ["{path[0]}/{1[0]}_1.fastq", "{path[0]}/{1[0]}_2.fastq"])
def collate_files(input_files, output_files):
    log.info("Collating paired fastq files: \n\t{} \n\t{}\n".format(input_files[0], input_files[1]))

@active_if(options.one_codex)
@active_if(options.paired)
@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}.assembled.fastq", extras)
def pear_fastq_files(input_files, output_file, extras):

    log.info("Starting pear run on %s and %s", input_files[0], input_files[1])

    output_file = re.sub(r"\.assembled\.fastq", "", output_file)
    args = ["pear", "-f", input_files[0], "-r", input_files[1], "-o", output_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("PEAR","Merging paried files with PEAR failed: \n{}".format(output))
        raise SystemExit

@active_if(options.one_codex)
@active_if(options.paired)
@transform(pear_fastq_files, suffix(".fastq"), ".fastq", extras)
def upload_paired_to_one_codex(input_file, output_file, extras):

    args = ["onecodex", "upload", input_file]
    log.info("uploading %s to One Codex", input_file)

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("One Codex","uploading peared files to One Codex failed \n{}".format(output))
        raise SystemExit

@active_if(options.one_codex)
@active_if(not options.paired)
@transform(input_files, suffix(".fastq"), ".fastq", extras)
def upload_to_one_codex(input_file, output_file, extras):

    args = ["onecodex", "upload", input_file]
    log.info("uploading %s to One Codex", input_file)

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("One Codex","uploading to One Codex failed \n{}".format(output))
        raise SystemExit

# paired alignment
@active_if(options.paired)
@transform(collate_files, formatter("([^/]+)_[12].fastq$"), options.output+"{1[0]}-tophat-results/accepted_hits.bam", options, extras)
def tophat_align_paired(input_files, output_file, options, extras):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    args = ["tophat2", "-G", options.gtf,"-p", str(options.cores), "-o", output, options.index, input_files[0], input_files[1]]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        extras.report_error("tophat_paired","tophat2 paired run failed:\n{}".format(output))
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# pre trimming
@active_if(options.trim)
@transform(input_files, suffix(".fastq"), ".trimmed.fastq", options, extras)
def trim_fastq_files(input_file, output_file, options, extras):

    # trim it
    args = "seqtk trimfq -e {} {} > {}".format(options.trim_val, input_file, output_file)
    p = subprocess.Popen(args, shell=True)
    p.wait()
    
    if p.returncode != 0: 
        log.warn("SeqTK failed trimming %s", input_file)
        extras.report_error("seqTK trimming", "failed")
        raise SystemExit

@active_if(options.trim)
@transform(trim_fastq_files, formatter(), options.output+"{basename[0]}-tophat-results/accepted_hits.bam", options, extras)
def tophat_align_trimmed(input_file, output_file, options, extras):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    log.info("Starting tophat2 run on trimmed %s", input_file)
    args = ["tophat2", "-G", options.gtf,"-p", options.cores, "-o", output, options.index, input_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn(output)
        extras.report_error("tophat_unpaird","tophat2 unpaired failed: \n{}".format(output))
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

    # delete the trimmed fastq file
    log.info("Deleting the trimmed file")
    os.unlink(input_file)

# unpaired alignment function
@active_if(not options.paired)
@active_if(not options.trim)
@transform(input_files, formatter(), options.output+"{basename[0]}-tophat-results/accepted_hits.bam", options, extras)
def tophat_align_unpaired(input_file, output_file, options, extras):

    # we cant to have tophat write results to output+filename
    output = os.path.dirname(output_file)
    log.info("Starting tophat2 run on %s", input_file)
    args = ["tophat2", "-G", options.gtf,"-p", options.cores, "-o", output, options.index, input_file]

    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn(output)
        extras.report_error("tophat_unpaird","tophat2 unpaired failed: \n{}".format(output))
        raise SystemExit

    # print output
    log.info("tophat2 output:")
    log.info(output)

# get rid of stupid accepted_hits.bam file
@transform([tophat_align_unpaired, tophat_align_paired, tophat_align_trimmed],
            formatter(r".*/([^/]+)-tophat-results/accepted_hits.bam$"), 
            ''.join([options.output, "{1[0]}", ".bam"]), extras )
def rename_accepted_hits(input_file, output_file, extras):

    try:
        os.rename(input_file, output_file)
    except OSError:
        extras.report_error("rename_accepted_hits","Renaming {} to {} failed".format(input_file, output_file),log)

# both the tophat functions can feed into here
@transform(rename_accepted_hits, suffix(".bam"),".sorted.bam", extras)
def sort_bam(input_file, output_file, extras):
    log.info("Sorting %s ", input_file)

    # hacky
    output_file = re.sub(r"\.bam", "", output_file)

    if subprocess.call(["samtools-rs", "rocksort", "-@", "8", "-m", "16G", input_file, output_file]):
        log.warn("bam sorting %s failed, exiting", input_file)
        extras.report_error("sort_bam","bam sorting {} failed".format(input_file))
        raise SystemExit
    
    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@active_if(options.bw)
@transform(sort_bam, suffix(".sorted.bam"), ".bed", options.output, extras)
def bam_to_bed(input_file, output_file, output, extras):

    log.info("Converting %s to a bed file", input_file)
    if subprocess.call("bamToBed -i {} > {}".format(input_file, output_file), shell=True):
        log.warn("bam to bed conversion of %s failed, exiting", input_file)
        extras.report_error("bam_to_bed","bam to bed conversion of {} failed".format(input_file))
        raise SystemExit

    # now we can move sorted bam to output
    file_name = os.path.basename(input_file)
    new_name = os.path.join(output, file_name)
    os.rename(input_file, new_name)

@active_if(options.bw)
@transform(bam_to_bed, suffix(".bed"), ".bg", options.size, extras)
def bed_to_bg(input_file, output_file,  size_file, extras):
    log.info("Converting %s to a genome coverage file", input_file)
    
    command = "genomeCoverageBed -bg -split -i {} -g {} > {}".format( input_file, size_file, output_file)
    if subprocess.call(command, shell=True):
        log.warn("bed to coverage conversion of %s failed, exiting", input_file)
        extras.report_error("bed_to_bg","bed to bg conversion of {} failed".format(input_file))
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

@active_if(options.bw)
@transform(bed_to_bg, suffix(".bg"), ".bw", genome, options.output)
def bg_to_bw(input_file, output_file, genome, output):
    log.info("Creating bigwig file from bg:  %s", input_file)
    command = "bedGraphToBigWig {} {} {}".format( input_file, size_file, output_filet)
    
    if subprocess.call(command, shell=True):
        log.warn("bg to bw conversion of %s failed, exiting", input_file)
        extras.report_error("bg_to_bw","bg to bw conversion of {} failed".format(input_file))
        raise SystemExit

    log.info("Deleting old file %s", input_file)
    os.unlink(input_file)

# run the pipeline
cmdline.run(options)
