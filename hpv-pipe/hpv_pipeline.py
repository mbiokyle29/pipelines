#!/usr/bin/env python
"""
hpv_pipeline.py
Kyle McChesney
"""

import logging
import operator
import os.path as path
import time
from collections import defaultdict as ddict
from ftplib import FTP

import ruffus.cmdline as cmdline
import sh
from ruffus import *

from hpv_utils import read_srx_list, read_srr_file

parser = cmdline.get_argparse(description='HPV pipeline')

parser.add_argument("--srx-file", "-s",
                    help="List file with desired SRX nums")

parser.add_argument("--srr-file", "-r",
                    help="List file with SRX: SRR,SRR,SRR info")

parser.add_argument("--index", "-i",
                    help="Path to BWA index /path/basename format",
                    required=True)

parser.add_argument("--fasta", "-f",
                    help="Path to fasta file (used to make index)",
                    required=True)

parser.add_argument("--cores", "-p",
                    help="Number of cores to run on", default=5)

parser.add_argument("--output", "-o", default="./",
                    help="Output / Working directory for everything")

# build ops, setup log
options = parser.parse_args()

### LOGGING
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

# file log
time_stamp = str(time.time()).replace(".","")
log_file = options.log_file if options.log_file else path.join(options.output,"{}.{}.{}".format("hpv_pipeline",time_stamp,"log"))
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

# set output to fp
output = path.abspath(options.output)

# set up the ftp connection
FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"
EMAIL = "kgmcchesney@wisc.edu"  # SRA ftp requires email as password
ftp_conn = FTP(FTP_ROOT)
ftp_conn.login(passwd=EMAIL)

# sh.py commands
try:
    fqd = sh.fastq_dump.bake("-O", options.output, "-A")
    bwa = sh.bwa.bake("mem", "-t", options.cores, options.index)
    stv = sh.samtools.bake("view", "-F", "4", "-bT", options.fasta)
    stvb = sh.samtools.bake("view", "-h")
    sts = sh.samtools.bake("sort")
    wcl = sh.wc.bake("-l")

except ImportError as e:
    log.warn("Import error occured: %s", e)
    log.warn("Please ensure: fastq-dump, bwa and samtools are installed and in your path")
    raise SystemExit

# get input data, either from srx or srr source
if options.srr_file:
    srr_data = read_srr_file(options.srr_file)

elif options.srx_file:
    srxs = build_srx_list(options.srx_file)
    srr_data = build_srr_list(srxs, ftp_conn)

else:
    log.error("Error: One of --srx-file or --srr-file must be supplied!")
    raise SystemExit

# determine the full list of SRR files we need
# accounting for paired end
total_srr_list = reduce(lambda x,y: x + y, [srr_data[z] for z in srr_data.keys()])
total_srr_list = map(lambda x: options.output + x +".fastq", total_srr_list)
assert len(set(total_srr_list)) == len(total_srr_list)

@originate(total_srr_list, fqd, output)
def fake_download_fastq_files(fastq_file, fqd, output):

    srr_file = fastq_file.replace(".fastq", ".sra")
    sh.touch(srr_file)

@transform(fake_download_fastq_files, suffix(".fastq"), ".sam")
def align_with_bwa(input_fastq, output_sam, fasta):

    srr_number = path.splitext(path.basename(srr_file))[0]
    log.info("Fastq dumping: %s (%s)", srr_file, srr_number)
    assert srr_number.startswith("SRR")
    fqd(srr_number)

    log.info("Aligning %s with %s", input_fastq, bwa)
    bwa_log = output_sam + ".log"
    bwa(input_fastq, _out=output_sam, _err=bwa_log)

@transform(align_with_bwa, suffix(".sam"), ".bam")
def sam_to_bam(sam_file, bam_file):
    
    log.info("Converting %s to a bam", sam_file)
    stv(sam_file, _out=bam_file)

@transform(sam_to_bam, suffix(".bam"), ".sorted.bam")
def bam_sort(bam_file, sorted_bam_file):

    log.info("Sorting %s", bam_file)
    sorted_bam_file = bam_file.replace(".sorted.bam", ".sorted")
    sts(bam_file, sorted_bam_file)

@transform(bam_sort, suffix(".sorted.bam"), ".sorted.sam")
def bam_to_sam(bam_file, sam_file):

    log.info("Converting %s to a sorted sam file", bam_file)
    stvb(bam_file, _out=sam_file)

@transform(bam_to_sam, suffix(".sorted.sam"), ".data.out")
def parse_sam_file(sam_file, data_file):

    scores = ddict(int)
    with open(sam_file, "r") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            genome = line.rstrip().split("\t")[2]
            scores[genome] += 1

    sorted_scores = sorted(scores.items(), reverse=True,
                           key=operator.itemgetter(1))
    file_root = sam_file.replace(".sorted.sam", "")
    
    fastq_file = file_root + ".fastq"
    bam_file = file_root + ".bam"
    sbam_file = file_root + ".sorted.bam"

    lines = wcl(fastq_file)
    num_lines = int(lines.split()[0])
    num_reads = num_lines / 4

    with open(data_file, "w+") as fh:
        fh.write("### DATA REPORT FOR {} ###\n".format())
        fh.write("fastq lines: \n{}\n".format(lines))
        fh.write("reads: {}\n".format(str(num_reads)))
        fh.write("\n### Genome Hit Data ###")
        for genome,score in sorted_scores:
            fh.write("{}\t{}\n".format(genome, str(score)))

    os.unlink(fastq_file)
    os.unlink(bam_file)
    os.unlink(sbam_file)

# run the pipelined
cmdline.run(options)
