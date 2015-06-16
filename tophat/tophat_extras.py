import os, glob, re, subprocess

# :)
MILLION = float(1000000)

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