#!/usr/bin/env python
"""
hpv_pipeline.py
Kyle McChesney
"""
 
import logging as log

import ruffus.cmdline as cmdline
import sh
from ruffus import *

from utils import make_fastq_list

parser = cmdline.get_argparse(description='HPV pipeline')
parser.add_argument("--multi-fasta", 
                    help="Mutli fasta file with all HPV genomes", 
                    required=True)

parser.add_argument("--dir", 
                    help="Fullpath to the directory where the FASTQ reads are located", 
                    required=True)

parser.add_argument("--cores", "-p",
                    help="Number of cores to run on", default=5)



options = parser.parse_args()
log.basicConfig(level=log.INFO)

fastq_files = make_fastq_list(options.dir)
log.info("The following fastq files were found:")
for fastq in fastq_files:
    log.info(fastq)

try:
    
    # get and configure rsem callers
    from sh import rsem_calculate_expression as rce

    rce = rce.bake("-p", str(options.cores), "--no-bam-output",
                   "--keep-intermediate-files")

    log.info("RSEM initalized:")
    log.info(rce)

except ImportError as e:
    log.warn("Is RSEM installed?")
    log.warn(e)
    raise SystemExit

@collate(fastq_files, formatter("([^/]+)_[12].fastq$"), ["{path[0]}/{1[0]}_1.fastq", "{path[0]}/{1[0]}_2.fastq"])
def collate_files(input_files, output_files):
    log.info("Collating paired fastq files: \n\t{} \n\t{}\n".format(input_files[0], input_files[1]))

@transform(collate_files, suffix(".fastq"), "", options, rce)
def rsem_calculate_expression(input_files, output_file, options, rce):

    output_file = output_file.replace("_1", "")
    log.info("Running rsem calculate expression on %s", input_files)
    log.info("CMD: %s %s %s %s %s", rce, input_files[0], input_files[1], options.hg_index, output_file)
    
    rce(input_files[0], options.hg_index, output_file)


# run the pipelined
cmdline.run(options)
