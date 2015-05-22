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

# parse the args
options = parser.parse_args()

#  set up logs
logger, logger_mutex = cmdline.setup_logging (__name__, options.log_file, options.verbose)


input_files = ["test1.fastq","test2.fastq"]
tophat_opts = {
	'cores': 8,
	'gtf': "/data/refs/hg19/hg19.gtf",
	'index': "/data/refs/hg19/hg19",
	'output': "/data/results/"
}

@transform(input_files, suffix(".fastq"),".bam", tophat_opts)
#@mkdir(input_files,)
def align_with_tophat(input_file, output_file, tophat_opts):
	output = tophat_opts.get('output')
	tophat_command = "tophat2 -p {0} --GTF {1} -o {2} {3} {4}".format( 
					tophat_opts.get('cores'),
					tophat_opts.get('gtf'),
					output + input_file, 
					tophat_opts.get('index'),
					input_file)
	print tophat_command





# run the pipeline
pipeline_printout_graph ( "flowchart.svg")
cmdline.run (options)
