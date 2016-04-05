#!/usr/bin/env python
"""
deseq_pipeline.py
Kyle McChesney

An RNAseq pipline for differential expression
- Alignment with RSEM
- DE with DESEQ
"""
import logging
import os
import subprocess
import time

import pandas as pd
import ruffus.cmdline as cmdline
import sh
from ruffus import *

from utils import guess_simple_design_matrix, log_line, make_fastq_list, send_report

parser = cmdline.get_argparse(description='RSEM and deseq2 pipeline')
parser.add_argument("--dir", 
                    help="Fullpath to the directory where the FASTQ reads are located", 
                    required=True)

parser.add_argument("--index",
                    help="Fullpath to the bowtie2 index in: /full/file/path/basename form",
                    default="/data/refs/hg19/hg19-rsem")

parser.add_argument("--output",
                    help="Fullpath to output directory",
                    default="./")

parser.add_argument("--name",
                    help="Optional experiment name",
                    dest="exp_name", default="rsem-deseq")

parser.add_argument("--to",
                    help="List of emails to send final reports too",
                    nargs="+", default=["mbio.kyle@gmail.com"])

parser.add_argument("--cores", "-p",
                    help="Number of cores to run on", default=5)

parser.add_argument("--size", 
                    help="Fullpath to size file",
                    default="/data/refs/hg19/hg19.sizes")

parser.add_argument("--no-bam",
                    help="Disable bam output on rsem",
                    default=False, action='store_true')

parser.add_argument("--no-de",
                    help="Optionally skip deseq2 step (Maybe really complex samples)",
                    default=False, action='store_true')

parser.add_argument("--no-wig",
                    help="Optionally skip wig generation",
                    default=False, action='store_true')

parser.add_argument("--design-file",
                    help="DESeq2 design matrix file (TSV) samples & treatments")

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

    if options.no_bam:
        rce = rce.bake("-p", str(options.cores), "--no-bam-output",
                       _out=log_line, _err_to_out=True)
    else:
        rce = rce.bake("-p", str(options.cores), "--output-genome-bam", 
                       _out=log_line, _err_to_out=True)
    
    rbw = rbw.bake(_out=log_line, _err_to_out=True)
    wtbw = wtbw.bake(_out=log_line, _err_to_out=True)

    log.info("RSEM initalized:")
    log.info(rce)
    log.info(rbw)
    log.info(rgdm)

except ImportError as e:
    log.warn("Is RSEM installed?")
    log.warn(e)
    raise SystemExit


##
## Globals :/
###
fastq_files = make_fastq_list(options.dir)
log.info("The following fastq files were found:")
for fastq in fastq_files:
    prog_log.info(fastq)

if options.design_file is None and not options.no_de:
    options.design_file = guess_simple_design_matrix(fastq_files, options.dir)

@transform(fastq_files, suffix(".fastq"), ".genes.results", options, rce)
def rsem_calculate_expression(input_file, output_file, options, rce):
    
    output_file = output_file.replace(".genes.results", "")

    log.info("Running rsem calculate expression on %s", input_file)
    log.info("CMD: %s %s %s %s", rce, input_file, options.index, output_file)
    
    rce(input_file, options.index, output_file)

@merge(rsem_calculate_expression, os.path.join(options.output, options.exp_name+".counts.xlsx"), options)
def parse_cnt_file(expression_files, output_file, options):

    # we are getting genome bams but we want the .cnt files in the .stat dir
    output_fp = os.path.abspath(options.dir)
    df_idx = [os.path.basename(x).split(".")[0] for x in expression_files]
    stat_dir_files = [os.path.join(output_fp, x+".stat")  for x in df_idx]
    cnt_files = [os.path.join(y, x+".cnt") for y,x in zip(stat_dir_files, df_idx)] 
    
    df_cols = ["Unalignable", "Alignable", "Multi-filtered", "Total Reads"]
    counts_df = pd.DataFrame(index=df_idx, columns=df_cols)

    for cnt_file in cnt_files:
        with open(cnt_file, "r") as fh:
            reads = fh.next().rstrip().split()
        row = {"Unalignable": reads[0], "Alignable": reads[1],
               "Multi-filtered": reads[2], "Total Reads": reads[3]}
        row_idx = os.path.basename(cnt_file).split(".")[0]
        counts_df.loc[row_idx] = pd.Series(row)

    counts_df.to_excel(output_file, "Alignment Stats")
    

@active_if(not options.no_wig and not options.no_bam)
@transform(rsem_calculate_expression, suffix(".genes.results"), ".wig", options, rbw, output_dir = options.output)
def rsem_bam_to_wig(input_file, output_file, options, rbw):
    
    input_file = input_file.replace(".genes.results",".bam")
    plot_name = output_file.replace(".wig", "")

    log.info("Running rsem bam 2 wig on %s", input_file)
    log.info("CMD: %s %s %s %s", rbw, input_file, output_file, plot_name)

    rbw(input_file, output_file, plot_name)

@active_if(not options.no_wig and not options.no_bam)
@transform(rsem_bam_to_wig, suffix(".wig"), ".bw", options, wtbw, output_dir = options.output)
def wig_to_big_wig(input_file, output_file, options, wtbw):

    log.info("Converting %s to a big wig", input_file)
    log.info("CMD: %s %s %s %s", wtbw, input_file, options.size, output_file)

    wtbw(input_file, options.size, output_file)

@merge(rsem_calculate_expression, os.path.join(options.output, "{}{}".format(options.exp_name, "-counts.mtx")), options, rgdm)
def rsem_generate_data_matrix(input_files, output_file, options, rgdm):
    
    ordered_samples = sorted(input_files)
    log.info("Running rsem generate data matrix")
    log.info("CMD: %s %s > %s", rgdm, " ".join(input_files), output_file)
    rgdm(*ordered_samples, _out=output_file)


@follows(rsem_generate_data_matrix, wig_to_big_wig, parse_cnt_file,
         mkdir(options.output+"bams/", options.output+"vis/",
               options.output+"counts/", options.output+"stat/"))
@merge([rsem_generate_data_matrix, parse_cnt_file], "out.email", options)
def clean_up(files, email_name, options):

    log.info("Preparing to clean up")
    mv = sh.mv
    move_tups = []

    stat_glob = options.dir+"*.stat/*"
    stat_dir_glob = options.dir+"*.stat/"
    stat_dest = options.output+"stat/"

    count_glob = options.dir+"*.genes.results"
    isoform_glob = options.dir+"*isoforms.results"
    count_dest = options.output+"counts/"

    log.info("Moving data files")
    
    if not options.no_bam:
        gbam_glob = options.dir+"*.genome.sorted.bam"
        gbam_dest = options.output+"bams/"
        move_tups.append((gbam_glob, gbam_dest))

    if not options.no_bam and not options.no_wig:
        wig_glob = options.output+"*.wig"
        bw_glob = options.output+"*.bw"
        vis_dest = options.output+"vis/"
        move_tups.append((wig_glob, vis_dest))
        move_tups.append((bw_glob, vis_dest))

    for glob, dest in move_tups:
        mv(sh.glob(glob), dest)

    mv(sh.glob(count_glob), count_dest)
    mv(sh.glob(isoform_glob), count_dest)
    mv(sh.glob(stat_glob), stat_dest)
    
    log.info("Deleting junk files")

    if not options.no_bam:
        bai_glob = options.dir+"*.bai"
        ubam_glob = options.dir+"*.genome.bam"
        transcript_glob = options.dir+"*transcript*"
        [sh.rm(sh.glob(g)) for g in [bai_glob, ubam_glob, transcript_glob]]

    sh.rmdir(sh.glob(stat_dir_glob))

    report_files = [files[0], files[1], log.handlers[0].baseFilename]
    subject = "RSEM/deseq2 pipeline"

    if options.exp_name != "rsem-deseq":
        subject += " - {}".format(options.exp_name)

    log.info("Sending report")
    send_report(list(options.to), subject, report_files)
    log.info("Run complete! Congradulations!")


# run the pipelined
cmdline.run(options)
