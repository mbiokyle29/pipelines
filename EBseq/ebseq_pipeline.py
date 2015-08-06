#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for all things tophat

"""

# ruffus imports
from ruffus import *
import ruffus.cmdline as cmdline
import logging, time, os, subprocess
from ebseq_extras import EbseqExtras 

# EMAIL
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

parser = cmdline.get_argparse(description="This is a pipeline for RSEM alignment and subsequent DE gene analysis with EBseq")

# Program arguments
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--cores", help="Number of cores to run RSEM on", default='10')
parser.add_argument("--index", help="Fullpath to the RSEM index in: /full/file/path/basename form", 
                    default="/data/refs/hg19/hg19-RSEM", required=True)

parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--conf", help="Fullpath to conf tsv file, (<sample><frag-mean><frag-sd><cond>)",
                    required=True)

parser.add_argument("--fdr", help=" false discovery rate to use for DE", default=0.05)

# reporting
parser.add_argument("--emails", help="Emails to send results too", default="mbio.kyle@gmail.com", nargs="+")

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
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("ebseq_pipeline",time_stamp,"log"))
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
log.info("Starting RSEM -> EBseq pipeline")

# build the extras
# generate the sample recs
extras = EbseqExtras(log)
extras.read_configuration(options.conf)

#
# Plan
# fastqs --> each { rsem-calc-express | rsem-plot-model }
#   \_ group by cond --> generate-data-matrix --> rsem-run-ebseq --> rsem-control-fdr
#   \_ res + graphs :: EMAIL

# get the input files together
input_files = extras.gen_fastq_list()
log.info("Samples are: ")
for file in input_files:
    log.info(file)

@transform(input_files, formatter(), options.output+"{basename[0]}.genes.results", options, extras)
def rsem_align(input_file, output_file, options, extras):
    
    mean_len = extras.get_mean_length(input_file)
    output_file = output_file.replace("\.genes\.results", "")

    log.info("Running rsem calc exp on %s", output_file)
    command = ["rsem-calculate-expression", "-p", str(options.cores), "--calc-ci", input_file, options.index, output_file]

    run_cmd(command)

@transform(rsem_align, suffix(".genes.results"), ".pdf", options)
def plot_model(input_file, output_file, options):

    sample_name = input_file.replace("\.genes\.results", "")
    command = ["rsem-plot-model", sample_name, output_file]

    run_cmd(command)

@merge(rsem_align, "gene_exp.mtx", extras)
def generate_exp_matrix(input_files, matrix, extras):

    sample_list = extras.gen_sample_list()
    command = ["rsem-generate-data-matrix", sample_list, ">", matrix]

    run_cmd(command)

@transform(generate_exp_matrix, suffix(".mtx"), ".diff", extras)
def run_ebseq(input_file, output_file, extras):

    cond_str = extras.gen_cond_str()
    command = ["rsem-run-ebseq", input_file, cond_str, output_file]

    run_cmd(command)

@transform(run_ebseq, suffix(".diff"), ".sigdiff", options)
def fdr_correct(input_file, output_file, options):

    command = ["rsem-control-fdr", input_file, str(options.fdr), output_file]

    run_cmd(command)

# not a ruffus command
# general command runner
def run_cmd(cmd):

    log.info("Running: %s", cmd)
    try:
        output = subprocess.check_output(cmd)
	log.info(output)
    # log the error, report it via email and exit
    except subprocess.CalledProcessError as e:
        log.error("Command failed with error: %s", e.message)
        report_error(cmd)
        raise SystemExit

# run the pipeline
cmdline.run(options)
