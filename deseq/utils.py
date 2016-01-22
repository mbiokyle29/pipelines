import os
import logging
import sh
import glob

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