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
def make_sam_list(directory):

    sams = []

    # make sure the dir exists
    if not os.path.isdir(directory):
        log.warn( "%s is not a valid dir, exiting", directory)
        raise SystemExit 

    directory = os.path.abspath(directory)
    log.info("Reading %s for sams", directory)

    blob = os.path.join(directory, "*.sam")
    sams.extend(glob.glob(blob))
    
    # make sure we got stuff
    if len(sams) == 0:
        log.warn("Sam list is empty, exiting")
        raise SystemExit

    return sams

def log_line(line):
    prog_log.info(line.rstrip())

def send_report(to, subj, files):
    
    # Create a text/plain message
    email_body = []
    email_body.append("RNAseq variant calling pipeline results")
    email_body.append("Please direct any questions to kgmcchesney@wisc.edu")

    # msg object results
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

def is_valid_file(file):
    fp = os.path.abspath(file)
    return os.path.isfile(file)
