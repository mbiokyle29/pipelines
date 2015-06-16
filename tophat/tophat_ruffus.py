#!/usr/bin/env python
"""
Kyle McChesney

"""
from ruffus import *
import ruffus.cmdline as cmdline
import subprocess

parser = cmdline.get_argparse(description='WHAT DOES THIS PIPELINE DO?')

# add args here
parser.add_argument("--cores", help="Number of cores to run on (for ruffus and tophat)")
parser.add_argument("--gtf", help="Fullpath to the GTF reference file", required=False)
parser.add_argument("--index", help="Fullpath to the bowtie2 index in: /full/file/path/basename form")
parser.add_argument("--output", help="Fullpath to output directory", default="./")

# parse the args
options = parser.parse_args()
print(options)
#  set up logs
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)
input_files = ["test1.fastq","test2.fastq"]



@transform(input_files, suffix(".fastq"),".bam", options)
#@mkdir(input_files,)
def align_with_tophat(input_file, output_file, tophat_opts):
	output = tophat_opts.output
	tophat_command = "tophat2 -p {0} --GTF {1} -o {2} {3} {4}".format( 
					tophat_opts.cores,
					tophat_opts.gtf,
					output + input_file, 
					tophat_opts.index,
					input_file)
	print tophat_command


# run the pipeline
cmdline.run (options)