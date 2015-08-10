import os.path
import re

# don't use slots since we only have a few of these guys
class _sampleRec():

    def __init__(self, name, mean, std, condition):
        self.name = name
        self.mean = int(mean)
        self.std = int(std)
        self.condition = int(condition)

class EbseqExtras():

    def __init__(self, log):
        self.log = log
        self.samples = []
        self.conditions = {}

    def read_configuration(self, conf):
        if os.path.isfile(conf):
            try:
                with open(conf, "r") as fh:
                    for line in fh:
                        self._build_rec(line)

            except IOError as e:
                self.log.error("IOError thrown trying to read %s conf file, perhap permissions?", conf)
                raise SystemExit
        else:
            self.log.error("It appears %s does not exist", conf)
            raise SystemExit

    def _build_rec(self, line):
        # <sample><frag-mean><frag-sd><cond>
        rec = _sampleRec(*line.split("\t"))
        self.samples.append(rec)

        if rec.condition in self.conditions:
            self.conditions[rec.condition].append(rec)
        else:
            self.conditions[rec.condition] = [rec]

    def gen_fastq_list(self):
        results = []
        for sample in self.samples:
            results.append(sample.name)
        return results

    def gen_sample_list(self):

        sample_str = ""
        
        for cond in sorted(self.conditions.keys()):
            for rec in self.conditions[cond]:
                name = re.sub(r"\.fastq", ".genes.results", rec.name)
                sample_str += name+" " 

        return sample_str.rstrip()

    def get_mean_length(self, file):
        
        base = os.path.splitext(file)[0]
        for sample in self.samples:

            sample_base = os.path.splitext(sample.name)[0

            ]
            if base == sample_base:
                return sample.mean

        # if it wasnt found
        raise SystemError

    def gen_cond_string(self):
        # if conditions has {1}:[2], {2}:[2], {3}:[2]
        # we want 2,2,2
        cond_str = ""

        for condition in sorted(self.conditions.keys()):
            cond_str += str(len(self.conditions[condition]))+","

        return cond_str.rstrip(",")