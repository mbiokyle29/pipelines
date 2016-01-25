import os
import logging
import sh
import glob
import re
from collections import defaultdict

# EMAIL
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email.Utils import COMMASPACE, formatdate
from email import Encoders


prog_log = logging.getLogger("__program__")
log = logging.getLogger("__main__")

# not a ruffus task
def make_fastq_list(directory):

    fastqs = []

    # make sure the dir exists
    if not os.path.isdir(directory):
        log.warn( "%s is not a valid dir, exiting", directory)
        raise SystemExit 

    directory = os.path.abspath(directory)
    log.info("Reading %s for fastqs", directory)

    # see if there are any compressed files
    gz_blob = os.path.join(directory, "*.fastq.gz")
    gzs = glob.glob(gz_blob)

    for gz in gzs:
        log.info("gunzipping %s", gz)   

        # 0 == all good, and bool false
        if sh.gunzip(gz):
            log.warn("gunzipping %s failed, exiting", gz)
            raise SystemExit

    # now glob the fastqs
    blob = os.path.join(directory, "*.fastq")
    fastqs.extend(glob.glob(blob))
    
    # make sure we got stuff
    if len(fastqs) == 0:
        log.warn("Fastq list is empty, exiting")
        raise SystemExit

    return fastqs

def log_line(line):
    prog_log.info(line.rstrip())

def send_report(to, subj, files):
    
    # Create a text/plain message
    email_body = []
    email_body.append("RSEM / deseq2 pipeline results")
    email_body.append("Please direct any questions to kgmcchesney@wisc.edu")

    # msg object
    msg = MIMEMultipart()

    # header stuff
    # no one else cares but me!
    root  = "root@alpha-helix.oncology.wisc.edu"

    msg['Subject'] = subj
    msg['From'] = root
    msg['To'] = COMMASPACE.join(to)
    msg.attach( MIMEText("\n".join(email_body)) )

    # attatch the files
    for f in files:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(f,"rb").read() )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(f))
        msg.attach(part)
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(root, to, msg.as_string())
    s.quit()

def guess_simple_design_matrix(fastq_list, path):

    # the hope is that the files look like
    # [STUFF]_rep#.fastq
    groups = defaultdict(list)
    pattern = re.compile(r"(.+)[-_.]rep\d\.fastq")

    for fq in fastq_list:
        fq = os.path.basename(fq)
        match = pattern.match(fq)
        if match is not None:
            groups[match.group(1)].append(fq)

    potential_matrix = ["\tcondition"]
    for group in sorted(groups):
        group_samples = sorted(groups[group])
        for samp in group_samples:
            potential_matrix.append("{}\t{}".format(samp, group))
    
    log.info("Here is a best guess matrix file\n")
    prog_log.info("\n".join(potential_matrix))
    prog_log.info("\nDoes this look ok?")
    answer = raw_input("(Y/N) -> ")

    if answer == "Y":
        matrix_file = os.path.join(path, "deseq-conditions.mtx")
        with open(matrix_file, "w") as fh:
            for line in potential_matrix:
                fh.write(line)
        return matrix_file
    else:
        log.warn("Guessed matrix does not work, and none supplied, exiting")
        raise SystemExit