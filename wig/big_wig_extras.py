import os, glob, re, subprocess


class BigWigExtras():
    
    def __init__(self,log):
        self.log = log
        self.MILLION = float(1000000)
        self.mr_hash = {}

    # not a ruffus task
    def make_fastq_list(self, directorys):
        
        fastqs = []
        for directory in directorys:

            # make sure the dir exists
            if not os.path.isdir(directory):
                self.log.warn( "%s is not a valid dir, exiting", directory)
                raise SystemExit 

            directory = os.path.abspath(directory)    
            # see if there are any compressed files
            gz_blob = os.path.join(directory, "*.fastq.gz")
            gzs = glob.glob(gz_blob)
        
            for gz in gzs:
                self.log.info("gunzipping %s", gz)   
        
                # 0 == all good, and bool false
                if subprocess.call(["gunzip",gz]):
                    self.log.warn("gunzipping %s failed, exiting", gz)
                    raise SystemExit
        
            # now glob the fastqs
            blob = os.path.join(directory, "*.fastq")
            fastqs.extend(glob.glob(blob))
        
        # make sure we got stuff
        if not fastqs:
            self.log.warn("Fastq list is empty, exiting")
            raise SystemExit

        return fastqs
    # not a ruffus command
    def record_bowtie_output(self, output, input_file, stats_file):
        # first make sure the output isnt messed up
        # should have 12 lines
        output = output.split("\n")
        output.pop()

        if not len(output) == 12:
            self.log.warn("Bowtie output is corrupted")
            raise SystemExit

        # things that we NEED from this record
        # <sample-name>\t<total-reads>\t<million reads>\t<reads mapped>
        # NOTE reads mapped is a sum of 1 time reads and >1 times

        # grab total reads
        # 4th index: NNNN reads; of these:

        total_reads = int(re.match(r"^(\d+)", output[4]).group(0))
        mill_reads  = total_reads / self.MILLION

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

        # update hash
        base = os.path.splitext(os.path.basename(input_file))[0]
        self.mr_hash[base] = mill_reads

        self.log.info("bowtie results for %s", input_file)
        self.log.info(entry)

    # not a ruffus command
    def make_mill_reads(self, file):
        
        if not os.path.isfile(file):
            self.log.warn( "%s is not a valid file, exiting", file)
            raise SystemExit

        with open(file, "r") as fh:
            for line in fh:
                # 0 is name, 2 is MRM
                arr = line.split("\t")
                name = os.path.splitext(os.path.basename(arr[0]))[0]
                mrm = arr[2]
                self.mr_hash[name] = float(mrm)

    def get_mill_reads(self, base):
        return self.mr_hash[base]

    # not a ruffus task
    def make_bam_list(self, directory):

        # make sure the dir exists
        if not os.path.isdir(directory):
            self.log.warn( "%s is not a valid dir, exiting", directory)
            raise SystemExit 

        blob = os.path.join(directory, "*.bam")
        bams = [os.path.abspath(bam) for bam in glob.glob(blob)]

        return bams