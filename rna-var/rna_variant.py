#!/usr/bin/env python
"""
rna_variant.py
Kyle McChesney

An RNAseq pipline for variant calling.
This basically implements the following: https://www.broadinstitute.org/gatk/guide/article?id=3891

A few notes:
- Steps 4/5 are not done. They are marked as optional by GATK right now, but should be implemented
- Step 7 (filtering) is also unimplemented, as this dependent upon project specification, however 
  some filter should be applied at some time!
- This accepts sam files, but expects them to be STAR alignments (GATK thinks thats best)
  if using none STAR alignments their is a MAPQ remapping option you will need to change
- There are areas for improvements in the add_or_replace_readgroups step.
  A .tsv file should be read-in mapping sam files to a more appropriate rg id, platform, etc
"""
import logging
import os
import subprocess
import time
from functools import partial

import ruffus.cmdline as cmdline
import sh
from ruffus import *


from utils import log_line, make_sam_list, send_report, is_valid_file

parser = cmdline.get_argparse(description='RNAseq variant calling pipeline')
parser.add_argument("--dir", 
                    help="Fullpath to the directory where the sam alignments are located", 
                    required=True)

parser.add_argument("--picard", required=True,
                    help="Fullpath to the picard tools .jar file")

parser.add_argument("--gatk", required=True,
                    help="Fullpath to the GATK .jar file")

parser.add_argument("--output",
                    help="Fullpath to output directory",
                    default="./")

parser.add_argument("--fasta",
                    help="Reference fasta, must also be an .fai present",
                    required=True)

parser.add_argument("--name",
                    help="Optional experiment name",
                    dest="exp_name", default="rsem-deseq")

parser.add_argument("--to",
                    help="List of emails to send final reports too",
                    nargs="+", default=["mbio.kyle@gmail.com"])

parser.add_argument("--cores", "-p",
                    help="Number of cores to run on", default=5)

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
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("rna_variant",time_stamp,"log"))

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

log.info("Starting RNAvariant pipeline run")

# gather deps
try:
    picard = None
    gatk = None
    java_base = sh.java
    if is_valid_file(options.picard):
        picard = java_base.bake("-jar", options.picard)
    else:
        raise SystemExit("Picard .jar path is invalid")
    if is_valid_file(options.gatk):
        gatk = java_base.bake("-jar", options.gatk)
    else:
        raise SystemExit("Gatk .jar path is invalid")
    assert is_valid_file(options.fasta)

except ImportError as e:
    log.error("Failure in gathering required deps, making sure Java is installed and jar paths are correct")

## Globals :/
sam_files = make_sam_list(options.dir)
log.info("The following sam files were found:")
for sam in sam_files:
    prog_log.info(sam)

@transform(sam_files, suffix(".sam"), ".rg.bam", options, picard)
def add_or_replace_readgroups(sam_file, output_file, options, picard):

    # set up the base
    base_command = picard.bake("AddOrReplaceReadGroups", "I={}".format(sam_file),
                               "O={}".format(output_file), "SO=coordinate")
    
    # get a fake index for read grouping
    # this should really be pulled in via a .tsv file
    rg_id = sam_files.index(sam_file)
    rg_id = partial(("{}=" + str(rg_id)).format)
    final_cmd = base_command.bake(rg_id("RGID"), rg_id("RGLB"), rg_id("RGPL"), rg_id("RGPU"), rg_id("RGSM"))
    log.info("Running AddorReplaceReadGroups: \n{}".format(str(final_cmd)))
    final_cmd()

@transform(add_or_replace_readgroups, suffix(".rg.bam"), ".dup.rg.bam", options, picard)
def mark_duplicates(rg_bam, dedupped_bam, options, picard):

    # set up the base
    metrics = dedupped_bam.replace(".dup.rg.bam", ".dedupe.metrics")
    dedupe_cmd = picard.bake("MarkDuplicates", "I={}".format(rg_bam), "O={}".format(dedupped_bam),
                             "CREATE_INDEX=true", "VALIDATION_STRINGENCY=SILENT", "M={}".format(metrics))
    log.info("Running MarkDuplicates: \n{}".format(str(dedupe_cmd)))
    dedupe_cmd()

@transform(mark_duplicates, suffix(".dup.rg.bam"), ".split.bam", options, gatk)
def split_n_cigar_reads(input_bam, output_bam, options, gatk):

    cmd = gatk.bake("-T", "SplitNCigarReads", "-R", options.fasta,
                    "-I", input_bam, "-o", output_bam, "-rf", "ReassignOneMappingQuality",
                    "-RMQF", 255, "-RMQT", 60, "-U", "ALLOW_N_CIGAR_READS")
    log.info("Running split reads: \n{}".format(str(cmd)))
    cmd()

@transform(split_n_cigar_reads, suffix(".split.bam"), ".vcf", options, gatk)
def run_haplotyper(input_bam, output_vcf, options, gatk):
    cmd = gatk.bake("-T", "HaplotypeCaller", "-R", options.fasta, "-I", input_bam,
                    "-dontUseSoftClippedBases", "-stand_call_conf", 20.0, "-stand_emit_conf", 20.0,
                    "-o", output_vcf)
    log.info("Running HaplotypeCaller: \{}".format(str(cmd)))
    cmd()

# run the pipelined
cmdline.run(options)
