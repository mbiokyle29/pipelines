import os, glob, re

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

# not a ruffus command
def record_bowtie_output(output,input_file, stats_file, log):
    # first make sure the output isnt messed up
    # should have 12 lines
    output = output.split("\n")
    output.pop()

    if not len(output) == 12:
        log.warn("Bowtie output is corrupted")
        raise SystemExit

    # things that we NEED from this record
    # <sample-name>\t<total-reads>\t<million reads>\t<reads mapped>
    # NOTE reads mapped is a sum of 1 time reads and >1 times

    # grab total reads
    # 4th index: NNNN reads; of these:

    total_reads = int(re.match(r"^(\d+)", output[4]).group(0))
    mill_reads  = total_reads / MILLION

    ## reads mapped
    ## sum the starting number in 7th and 8th
    one_align = int(re.match(r"^\s+(\d+)", output[7]).group(0))
    mult_align = int(re.match(r"^\s+(\d+)", output[8]).group(0))
    total_aligned = one_align + mult_align

    # convert input_file name into a stats file, right the record and return MRM
    # so it can be stored in main
    # TODO this is nasty
    entry = "{}\t{}\t{}\t{}\n".format(input_file, total_reads, mill_reads, total_aligned)

    f = open(stats_file,"a") #opens file with name of "test.txt"
    f.write(entry)
    f.close() 

    log.info("bowtie results for %s", input_file)
    log.info(entry)

    # return MRM
    return mill_reads