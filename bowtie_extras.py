import os, glob

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

# not a ruffus command
def make_fastq_list(directory, log):
    
    # make sure the dir exists
    if not os.path.isdir(directory):
        log.warn( "%s is not a valid dir, exiting", directory)
        raise SystemExit 

    # glob the fastqs
    # this is kinda hacky right now (need re for optional .gz)
    blob = os.path.join(directory, "*.fastq*")
    fastqs = glob.glob(blob)

    # make sure we got stuff
    if not fastqs:
        log.warn("Fastq list is empty, exiting")
        raise SystemExit

    return fastqs