#!/usr/bin/env python
"""
hpv_pipeline.py
Kyle McChesney
"""

import logging as log
from ftplib import FTP

import ruffus.cmdline as cmdline
import sh
from ruffus import *

from hpv_utils import build_srx_list

FTP_ROOT = "ftp-trace.ncbi.nlm.nih.gov"
PATH_ROOT = "/sra/sra-instant/reads/ByExp/sra/SRX/SRX{}/SRX{}/"
EMAIL = "kgmcchesney@wisc.edu"  # SRA ftp requires email as password

parser = cmdline.get_argparse(description='HPV pipeline')
parser.add_argument("--reference-fasta", "-r" 
                    help="Mutli fasta file with all HPV genomes", 
                    required=True)

parser.add_argument("--srx-file", "-s" 
                    help="List file with desired SRX nums", 
                    required=True)

parser.add_argument("--cores", "-p",
                    help="Number of cores to run on", default=5)

parser.add_argument("--db", "-d",
                    help="SRAdb sqlite3 file",
                    required=True)

options = parser.parse_args()
log.basicConfig(level=log.INFO)

# set up sh.py
touch = sh.touch

# set up the ftp connection
ftp_conn = FTP(FTP_ROOT)
ftp_conn.login(passwd=EMAIL)

# read in the SRX ids
srx_ids = build_srx_list(opts.srx_file)

# run the pipelined
cmdline.run(options)
