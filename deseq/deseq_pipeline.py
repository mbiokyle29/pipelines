#!/usr/bin/env python
"""
deseq_pipeline.py
Kyle McChesney

An RNAseq pipline for differential expression
- Alignment with RSEM
- DE with DESEQ
"""

from ruffus import *
import ruffus.cmdline as cmdline
import subprocess
import logging
import sh
import time
import os
from rpy2.robjects.packages import importr
import rpy2.robjects as ro


# UTILS
from utils import make_fastq_list, log_line, send_report, guess_simple_design_matrix

parser = cmdline.get_argparse(description='RSEM and deseq2 pipeline')
parser.add_argument("--dir", help="Fullpath to the directory where the FASTQ reads are located", required=True)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form", default="/data/refs/AKATA-GFP/AKATA_GFP_RSEM")
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--name", help="Optional experiment name", dest="exp_name", default="rsem-deseq")
parser.add_argument("--to", help="List of emails to send final reports too", nargs="+", default=["mbio.kyle@gmail.com"])
parser.add_argument("--cores", "-p", help="Number of cores to run on", default=5)
parser.add_argument("--size", help="Fullpath to size file", default="/data/refs/AKATA-GFP/AKATA_GFP.size")

parser.add_argument("--no-de", help="Optionally skip deseq2 step (Maybe really complex samples)", default=False)
parser.add_argument("--no-wig", help="Optionally skip wig generation", default=False)
parser.add_argument("--design-file", help="DESeq2 design matrix file (TSV) samples & treatments")

options = parser.parse_args()

# Kenny loggins
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

prog_log = logging.getLogger("__program__")
prog_log.setLevel(logging.INFO)
prog_formatter = logging.Formatter('%(message)s')

# file log
time_stamp = str(time.time()).replace(".","")
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("deseq_pipeline",time_stamp,"log"))

file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(log_formatter)

prog_file_handler = logging.FileHandler(log_file)
prog_file_handler.setLevel(logging.INFO)
prog_file_handler.setFormatter(prog_formatter)

# console log
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(log_formatter)

prog_stream_handler = logging.StreamHandler()
prog_stream_handler.setLevel(logging.INFO)
prog_stream_handler.setFormatter(prog_formatter)

# set it all up
log.addHandler(file_handler)
log.addHandler(stream_handler)

prog_log.addHandler(prog_file_handler)
prog_log.addHandler(prog_stream_handler)

log.info("Starting DESEQ pipeline run")

try:
    
    # get and configure rsem callers
    from sh import rsem_calculate_expression as rce
    from sh import rsem_bam2wig as rbw
    from sh import rsem_generate_data_matrix as rgdm

    # get kent tools wig --> big wig
    from sh import wigToBigWig as wtbw

    rce = rce.bake("-p", str(options.cores), "--output-genome-bam", "--bowtie2", _out=log_line, _err_to_out=True)
    rbw = rbw.bake(_out=log_line, _err_to_out=True)
    wtbw = wtbw.bake(_out=log_line, _err_to_out=True)

    from sh import Rscript as rscript

    log.info("RSEM initalized:")
    log.info(rce)
    log.info(rbw)
    log.info(rgdm)

except ImportError as e:
    log.warn("Is RSEM installed?")
    log.warn(e)
    raise SystemExit

fastq_files = make_fastq_list(options.dir)

log.info("The following fastq files were found:")
for fastq in fastq_files:
    prog_log.info(fastq)

if options.design_file is None and not options.no_de:
    options.design_file = guess_simple_design_matrix(fastq_files, options.dir)

@transform(fastq_files, suffix(".fastq"), ".genome.sorted.bam", options, rce, output_dir = options.output)
def rsem_calculate_expression(input_file, output_file, options, rce):
    
    output_file = output_file.replace(".genome.sorted.bam", "")

    log.info("Running rsem calculate expression on %s", input_file)
    log.info("CMD: %s %s %s %s", rce, input_file, options.index, output_file)
    
    rce(input_file, options.index, output_file)

@active_if(not options.no_wig)
@transform(rsem_calculate_expression, suffix(".genome.sorted.bam"), ".wig", options, rbw, output_dir = options.output)
def rsem_bam_to_wig(input_file, output_file, options, rbw):
    plot_name = output_file.replace(".wig", "")

    log.info("Running rsem bam 2 wig on %s", input_file)
    log.info("CMD: %s %s %s %s", rbw, input_file, output_file, plot_name)

    rbw(input_file, output_file, plot_name)

@active_if(not options.no_wig)
@transform(rsem_bam_to_wig, suffix(".wig"), ".bw", options, wtbw, output_dir = options.output)
def wig_to_big_wig(input_file, output_file, options, wtbw):

    log.info("Converting %s to a big wig", input_file)
    log.info("CMD: %s %s %s %s", wtbw, input_file, options.size, output_file)

    wtbw(input_file, options.size, output_file)


@transform(rsem_calculate_expression, suffix(".genome.sorted.bam"), ".genes.results", output_dir = options.output)
def collect_gene_results(input_file, output_file):
    log.info("collecting count results for %s", input_file)

@merge(collect_gene_results, os.path.join(options.output, "{}{}".format(options.exp_name, "-counts.mtx")), options, rgdm)
def rsem_generate_data_matrix(input_files, output_file, options, rgdm):
    
    ordered_samples = sorted(input_files)
    log.info("Running rsem generate data matrix")
    log.info("CMD: %s %s > %s", rgdm, " ".join(input_files), output_file)
    rgdm(*ordered_samples, _out=output_file)

@active_if(not options.no_de)
@split(rsem_generate_data_matrix, "*.csv", options)
def run_deseq(matrix_file, csv_files, options):
    print "OK"

@follows(rsem_generate_data_matrix, wig_to_big_wig, 
         mkdir(options.output+"bams/", options.output+"vis/", options.output+"etc/",
               options.output+"counts/"))
@transform(rsem_generate_data_matrix, suffix(".mtx"), ".email", options)
def clean_up(matrix_file, email_name, options):

    log.info("Preparing to clean up")
    mv = sh.mv

    count_glob = options.output+"*.genes.results"
    count_dest = options.output+"counts/"

    gbam_glob = options.output+"*.genome.sorted.bam"
    gbam_dest = options.output+"bams/"
    
    wig_glob = options.output+"*.wig"
    bw_glob = options.output+"*.bw"
    vis_dest = options.output+"vis/"

    isoform_glob = options.output+"*isoforms.results"
    transcript_glob = options.output+"*transcript*"
    etc_dest = options.output+"etc/"

    bai_glob = options.output+"*.bai"
    ubam_glob = options.output+"*.genome.bam"

    log.info("Moving data files")
    mv(sh.glob(gbam_glob), gbam_dest)
    mv(sh.glob(count_glob), count_dest)
    mv(sh.glob(wig_glob), vis_dest)
    mv(sh.glob(bw_glob), vis_dest)
    mv(sh.glob(isoform_glob), etc_dest)
    
    log.info("Deleting junk files")
    [sh.rm(sh.glob(str(g))) for g in [bai_glob, ubam_glob, transcript_glob]]

    report_files = [matrix_file, log.handlers[0].baseFilename]
    subject = "RSEM/deseq2 pipeline"

    if options.exp_name != "rsem-deseq":
        subject += " - {}".format(options.exp_name)

    log.info("Sending report")
    send_report(list(options.to), subject, report_files)
    log.info("Run complete! Congradulations!")


# run the pipelined
cmdline.run(options)