import os, glob, re, subprocess, json, time
from jsonschema import validate
import smtplib
from email.mime.text import MIMEText


# not a ruffus command
def check_default_args(cores,index,output, log):
    
    def_opts = {
        "cores":10,
        "index":"/data/refs/hg19/hg19", 
        "output":"./"
    }

    running_opts = {
        "cores":cores,
        "index":index, 
        "output":output
    }

    for key in def_opts.keys():
        if def_opts[key] == running_opts[key]:
            log.warn(" {} parameter was not given using default: {}".format(key, def_opts[key]))

# not a ruffus task
def make_fastq_list(directory, log):
    
    # make sure the dir exists
    if not os.path.isdir(directory):
        log.warn( "%s is not a valid dir, exiting", directory)
        raise SystemExit 

    # see if there are any compressed files
    gz_blob = os.path.join(directory, "*.fastq.gz")
    gzs = glob.glob(gz_blob)

    for gz in gzs:
        log.info("gunzipping %s", gz)   

        # 0 == all good, and bool false
        if subprocess.call(["gunzip",gz]):
            log.warn("gunzipping %s failed, exiting", gz)
            raise SystemExit

    # now glob the fastqs
    blob = os.path.join(directory, "*.fastq")
    fastqs = [os.path.abspath(fastq) for fastq in glob.glob(blob)]

    # make sure we got stuff
    if not fastqs:
        log.warn("Fastq list is empty, exiting")
        raise SystemExit

    return fastqs

# not a ruffus task
def process_de_conf(file, log):

    # make sure the file is there
    if not os.path.isfile(file):
        log.warn( "%s is not a valid file, exiting", file)
        raise SystemExit 

    # schema to validate against
    de_conf_schema = {
        "type": "object",
        "properties": {
            "negative-condition": {
                "type":"object",
                "properties": {
                    "label": {"type": "string"},
                    "files": {"type": "array"}
                }
            },
            "positve-condition": {
                "type":"object",
                "properties": {
                    "label": {"type": "string"},
                    "files": {"type": "array"}
                }
            }
        }
    }
    with open(file) as json_data:
        conf = json.load(json_data)
        
        # make sure it matches
        try:
            validate(conf, de_conf_schema)

        except ValidationError:
            log.warn("%s is not valid DE conf json", conf)
            raise SystemExit

        # if its good return it
        return conf


# not a ruffus task
# send email
# the plan is to report the specific error
# and also append the the contents of the log for good measure
def report_error(task, message, log):


    # Create a text/plain message
    email_body = []
    email_body.append("Hello, Kyle\n")
    email_body.append("{} task failed with the following error: ".format(task))
    email_body.append(message)

    # grab the log file name from the log
    # we add the file handler first
    # so its here
    log_file = log.handlers[0].baseFilename
    email_body.append("\n#######################################################")
    email_body.append("#                PIPELINE LOG                         #")
    email_body.append("#######################################################")
    with open(log_file, "r") as log:
        for line in log:
            email_body.append(line.rstrip())

    msg = MIMEText("\n".join(email_body))

    # header stuff
    # no one else cares but me!
    root  = "root@alpha-helix.oncology.wisc.edu"
    me = "kgmcchesney@wisc.edu"
    subject = "Tophat pipeline failure report: {}".format(time.strftime("%d/%m/%Y"))
    
    msg['Subject'] = subject
    msg['From'] = root
    msg['To'] = me
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(root, [me], msg.as_string())
    s.quit()



def annotate_gene(gene):
    return "ANOOTATED"+gene