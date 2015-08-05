import os.path

# don't use slots since we only have a few of these guys
class _sampleRec():

    def __init__(self, name, mean, std, condition):
        self.name = name
        self.mean = mean
        self.std = std
        self.condition = condition

class EbseqExtras():

    def __init__(self, log):
        self.log = log
        self.samples = []

    def read_configuration(self, conf):
        if os.path.isfile(conf):
            try:
                with open(conf, "r") as fh:
                    for line in fh:
                        self._build_rec(line)

            except IOError as e:
                log.error("IOError thrown trying to read %s conf file, perhap permissions?", conf)
                # TODO return error maybe?
        else:
            log.error("It appears %s does not exist", conf)

    def _build_rec(self, line):
        # <sample><frag-mean><frag-sd><cond>
        sample, frag-mean, frag-sd, cond = line.split("\t")
        rec = _sampleRec(sample, frag-mean, frag-sd, cond)
        self.samples.append(rec)

    def gen_sample_list(self):

        # get each condition
        condition_indicies = {}
        for sample in self.samples:
            condition_indicies[sample.condition] = 1

        condition_str = {(index, array) for (index, []) in condition_indicies.keys() }

        for sample in self.samples:
            condition_str[sample.condition].append(sample.name+",")

        final_str = ""
        for condition_set in condition_str.keys():
            final_str.append([sample for sample in condition_str[condition_set]])

        return final_str