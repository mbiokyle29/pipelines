#!/usr/bin/env python
"""
Kyle McChesney

Ruffus pipeline for differential expression runs
Assuming we made bams earlier

"""

# ruffus imports
from ruffus import *
import ruffus.cmdline as cmdline

# custom functions
from tophat_extras import TophatExtras

# system imports
import subprocess, logging, os, re, time

# :) so i never have to touch excel
import pandas as pd

# for cummerbund
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

# EMAIL
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders

parser = cmdline.get_argparse(description='This pipeline provides a number of funtionalities for working with RNAseq data')

# Program arguments
parser.add_argument("--dir", help="Fullpath to the directory where the bams are located", required=True)
parser.add_argument("--cores", help="Number of cores to run cuffdiff on", default='10')
parser.add_argument("--output", help="Fullpath to output directory", default="./")
parser.add_argument("--size", help="Fullpath to size file")
parser.add_argument("--gtf", help="Fullpath to gtf file", required=True)
parser.add_argument("--de-conf", help="fullpath to differential expresssion configuration file", required=True)

parser.add_argument("--annotation-db", help="fullpath to the sqlite db file, <id><name><desc>")
parser.add_argument("--annotation-file", help="fullpath to a tsv file of gene annotations, will create sqlite db")

# reporting
parser.add_argument("--emails", help="Emails to send DE results too", default="mbio.kyle@gmail.com", nargs="+")

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
log_file = options.log_file if options.log_file else os.path.join(options.output,"{}.{}.{}".format("cuffdiff_pipeline",time_stamp,"log"))
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
log.info("Starting Cuffiff Run")

# pre checking
# create extas instance
extras = TophatExtras(log)

# grab the conf
conf = extras.process_de_conf(options.de_conf)
neg_files = []
pos_files = []


# checking the annotation situation :)
# we either got a db, got a db file and need to make a db, or do nothing
if options.annotation_file:
    extras.make_annotation_db(options.annotation_file, time_stamp, options.output)

elif options.annotation_db:
    extras.init_db(options.annotation_db)

input_files = extras.make_bam_list(options.dir)


@merge(input_files, os.path.join(options.output,"cuffdiff/gene_exp.diff",), options, extras)
def run_cuffdiff(input_files, output_file, options, extras):


    # make the output file
    output_dir = os.path.join(options.output,"cuffdiff/")

    # get the label
    neg_label = conf['negative-condition']['label']
    pos_label = conf['positive-condition']['label']
    label_string = "{},{}".format(neg_label, pos_label)

    # match the files in conf to the bams we have
    for conf_file in conf['negative-condition']['files']:

        conf_basename = os.path.splitext(conf_file)[0]
        # check the input files
        for bam_file in input_files:
            bam_basename = os.path.basename(bam_file).split(os.extsep)[0]
            
            # see if they match
            if conf_basename == bam_basename:
                log.info("%s matched %s in negative-condition", conf_file, bam_file)
                neg_files.append(bam_file)

                # remove it for optimization :) I think
                input_files.remove(bam_file)

     # match the files in conf to the bams we have
    for conf_file in conf['positive-condition']['files']:
        conf_basename = os.path.splitext(conf_file)[0]

        # check the input files
        for bam_file in input_files:
            bam_basename = os.path.basename(bam_file).split(os.extsep)[0]

            # see if they match
            if conf_basename == bam_basename:
                log.info("%s matched %s in positive-condition", conf_file, bam_file)
                pos_files.append(bam_file)

                # remove it for optimization :) I think
                input_files.remove(bam_file)


    # make sure both groups got something and that the inputfiles is emtpy
    if len(neg_files) == 0 or len(pos_files) == 0:
        log.warn("Cuffdiff file specification error!")
        extras.report_error("cuffdiff", "Cuffdiff file specification error!")
        raise SystemExit

    # call it
    log.info("Starting cuffdiff run")
    args = ["cuffdiff", "--FDR", "0.05", "-o", output_dir, "-L", label_string, "-p", options.cores, options.gtf, ','.join(neg_files), ','.join(pos_files)]
    try:
        output = subprocess.check_output(args, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        log.warn("Cuffdiff failed: %s", output)
        extras.report_error("cuffdiff","Cuffdiff failed: \n{}".format(output))
        raise SystemExit

    # print output
    log.info("cuffdiff output:")
    log.info(output)

@transform(run_cuffdiff, suffix(".diff"), ".xlsx", options, extras)
def write_excel_sheet(input_file, output_file, options, extras):
    
    # best message ever
    log.info("Writing gene expression results to spreadsheets")
    keep_cols = ["test_id", "locus", "sample_1", "sample_2", "status", "value_1",
                "value_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"]

    # read it in
    df = pd.read_table(input_file, sep="\t", usecols=keep_cols)

    # change the na,e of test_id
    keep_cols[0] = "gene_name"
    df.columns = keep_cols

    # annotate the genes
    if(options.annotation_db or options.annotation_file):
        df['gene_annotation'] = df['gene_name'].apply(extras.annotate_gene)

    # sort sig to the top
    df.sort(columns="significant", ascending=False, inplace=True, kind="heapsort")

    # write the sucker to a file
    writer = pd.ExcelWriter(output_file)
    df.to_excel(writer, sheet_name="Sheet1", index=False)
    writer.save()

    log.info("results written to: %s", output_file)

@transform(run_cuffdiff, suffix(".diff"), ".pdf", options)
def cummeRbund(input_file, output_file, options):

    # import the grapher and cummerBund
    grdevices = importr('grDevices')
    r_bund = importr("cummeRbund")
    r_plot = robjects.r('plot')


    # read in the diff results
    cuff = r_bund.readCufflinks(options.output)
    grdevices.pdf(file=input_file, width=10, height=10)

    r_plot(r_bund.dispersionPlot(r_bund.genes(cuff)))
    r_plot(r_bund.csBoxplot(r_bund.genes(cuff),replicates=True))
    r_plot(r_bund.csDendro(r_bund.genes(cuff),replicates=True))
    r_plot(r_bund.csBoxplot(r_bund.genes(cuff)))
    r_plot(r_bund.csDistHeat(r_bund.genes(cuff)))
    r_plot(r_bund.csDistHeat(r_bund.genes(cuff), replicates=True))
    r_plot(r_bund.PCAplot(r_bund.genes(cuff),"PC1","PC2"))
    r_plot(r_bund.PCAplot(r_bund.genes(cuff),"PC1","PC2",replicates=True))
    r_plot(r_bund.PCAplot(r_bund.genes(cuff),"PC2","PC2"))
    r_plot(r_bund.PCAplot(r_bund.genes(cuff),"PC3","PC2",replicates=True))

    # close the dev
    grdevices.dev_off()

@merge([write_excel_sheet, cummeRbund], ".email", options, input_files, extras)
def report_success(input_files, output_file, options, inputfiles, extras):
    
    log.info("Sending email report")
    
    # Create a text/plain message
    email_body = []
    email_body.append("Differential expression pipeline results:\n")
    email_body.append("The following bam files were compared:")
    
    email_body.append("negative condition:")
    for file in neg_files:
        email_body.append("- {}".format(file))

    email_body.append("positive condition:")
    for file in pos_files:
        email_body.append("- {}".format(file))

    email_body.append("\nThe results (xlsx spreadsheet) and pipeline log are attatched")
    email_body.append("Please direct any questions to kgmcchesney@wisc.edu")

    # msg object
    msg = MIMEMultipart()

    # header stuff
    # no one else cares but me!
    root  = "root@alpha-helix.oncology.wisc.edu"
    subject = "Cuffdiff Report for {} vs {} : {}".format(conf['negative-condition']['label'], conf['positive-condition']['label'], time.strftime("%d/%m/%Y"))

    msg['Subject'] = subject
    msg['From'] = root
    msg['To'] = COMMASPACE.join(options.emails)
    msg.attach( MIMEText("\n".join(email_body)) )

    # attatch the files
    for file in [input_files, log.handlers[0].baseFilename]:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(file,"rb").read() )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file))
        msg.attach(part)
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(root, options.emails, msg.as_string())
    s.quit()

# run the pipeline
cmdline.run(options)