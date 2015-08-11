#!/usr/bin/env python
"""
seq_pipe.py
Kyle McChesney

This script is the entry point to several different pipelines 
all of which have to do with NGS analysis.

Currently:

+ bowtie-alignment
+ tophat-alignment
+ cuffdiff
+ bigWig generation

"""

# ruffus imports
from ruffus import *
import ruffus.cmdline as cmdline

# custom functions
#from tophat_extras import TophatExtras

# system imports
import subprocess, logging, os, re, time
import pandas as pd
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

parser = cmdline.get_argparse(description='seq_pipe: A pipeline for performing various NGS analysis tasks')

# this is the only 'required' argument
# will control which pipeline is run
parser.add_argument("--analysis", help=" What type of analysis to perform",
                    choices = ["bowtie","tophat","cuffdiff","DE","bigWig"],
                    required=True)

# Program arguments
parser.add_argument("--input-dir", help="Fullpath to the directory where the input files are located", required=True)
parser.add_argument("--cores", help="Number of cores to run multi-threaded programs on", default='3')
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/hg19/hg19")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file")
parser.add_argument("--gtf", help="Fullpath to gtf file")
parser.add_argument("--paired", help="Indicates whether the reads in --dir are paired_end. MUST FOLLOW _1 _2 convention", default=False)

# optional arguments to control turning on and off tasks
parser.add_argument("--one-codex", help="Whether or not to upload each sample to one codex for metagenomic analysis", default=False)
parser.add_argument("--de-conf", help="fullpath to differential expresssion configuration file")
parser.add_argument("--annotation-db", help="fullpath to the sqlite db file, <id><name><desc>")
parser.add_argument("--annotation-file", help="fullpath to a tsv file of gene annotations, will create sqlite db")

# reporting
parser.add_argument("--emails", help="Emails to send DE results too", default="mbio.kyle@gmail.com", nargs="+")

# parse the args
options = parser.parse_args()

# Kenny loggins
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

# file log
time_stamp = str(time.time()).replace(".","")
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("seq_pipe",time_stamp,"log"))
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

# header logged
log.info("Seq Pipe Analysis starting")
log.info("Starting %s Analysis", options.analysis)
log.info("#####################################################################")

input_files = ["/home/kyle/Dev/pipelines/test-data/small.fastq"]
from tasks import upload_to_one_codex

# run the pipelined
cmdline.run(options)