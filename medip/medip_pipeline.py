#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for simple bowtie alignment
Does not support paired end reads

"""
import logging
import argparse
from rpy2.robjects import r
from rpy2.robjects.packages import importr

parser = cmdline.get_argparse(description='Perform the MEDIP analysis on a set of bam files')

# really just need bams i think
#parser.add_argument("--output", help="Fullpath to output directory", default="./")
#parser.add_argument("--cond-one-input", help=" Space seperated list of the condition one input files", required=True, nargs="+")
#parser.add_argument("--cond-one-treated", help=" Space seperated list of the condition one treated files", required=True, nargs="+")
#parser.add_argument("--cond-two-input", help=" Space seperated list of the condition two input files", required=True, nargs="+")
#parser.add_argument("--cond-two-treated", help=" Space seperated list of the condition two treated files", required=True, nargs="+")

# parse the args
options = parser.parse_args()

# logging
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO,)
log = logging.getLogger(__name__)
log.info("Starting MEDIPs run")


def roi_to_bed(roi_file):
    # going for a pretty simple bed here
    # chrom start stop name score
    # score will just follow the difference since they are s
    basename = os.path.basename(roi_file)
    bedfile  = ''.join(basename, ".bed")
    with open(roi_file, "r") as roi:
        with open(bedfile, "a+") as bed:
            header = next(roi)
            for line in roi:
                [chr, start, stop, first_mean, second_mean] = line.split("\t")
                score = (second_mean - first_mean) if second_mean > first_mean else (first_mean - second_mean)
                bed.write("{}\t{}\t{}\t{}\t{}".format(chr,start, stop, name, score))
