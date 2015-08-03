import os, glob, re, subprocess
import smtplib
from email.mime.text import MIMEText

class BowtieExtras():
    
    def __init__(self,log):
        self.log = log
        self.MILLION = float(1000000)
        self.mr_hash = {}

    # not a ruffus command
    def check_default_args(self, cores,index,output):
        
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
                self.log.warn(" {} parameter was not given using default: {}".format(key, def_opts[key]))

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
    def record_bowtie_output(self, output, input_file, stats_file, bowtie):
        
        output = output.split("\n")

        # set here so we can assign in block
        entry = None
        total_reads = 0
        mill_reads = 0 
        total_aligned = 0

        if bowtie:
            # lines 3 4 and 5 are what we want
            # they have a # at the start which is nice i guess

            if not len(output) == 9:
                self.log.warn("Bowtie output is corrupted")
                raise SystemExit

            total_reads = int(re.match(r"#\sreads\sprocessed:\s(\d+)$",output[3]).group(0))
            mill_reads = total_reads / self.MILLION
            total_aligned = int(re.match(r"#\sreads\swith\sat\sleast\sone\sreported\salignment:\s(\d+)\s", 
                                         output[4]).group(0))

        else:
            # first make sure the output isnt messed up
            # should have 12 lines
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

    def get_mill_reads(self, base):
        return self.mr_hash[base]

    # not a ruffus task
    # send email
    # the plan is to report the specific error
    # and also append the the contents of the log for good measure
    def report_error(self, task, message):

        # Create a text/plain message
        email_body = []
        email_body.append("Hello, Kyle\n")
        email_body.append("{} task failed with the following error: ".format(task))
        email_body.append(message)

        # grab the log file name from the log
        # we add the file handler first
        # so its here
        log_file = self.log.handlers[0].baseFilename
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
