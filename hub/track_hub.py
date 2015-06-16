import os, subprocess
class TrackHub:

    # init
    def __init__(self,root_dir, genome):
        
        # settin
        self.root_dir = root_dir
        self.genome = genome
        self.genome_dir = os.path.join(self.root_dir,self.genome)
        
        # make sure root exists, otherwise make
        if not os.path.isdir(self.root_dir):
            os.mkdir(self.root_dir)

        # make the subdirectory based on genome
        os.mkdir(self.genome_dir)

    def make_2bit_file(self, fasta_file):

        # make sure fasta exists
        if not os.path.isfile(fasta_file):
            print "Invalid fasta file"
            raise SystemExit

        # do it
        bit_file = os.path.join(self.genome_dir, self.genome, ".2bit")
        bit_command = ["faToTwoBit", fasta_file, bit_file]
        if subprocess.call(bit_command):
            print "Making 2bit file failed"
            raise SystemExit

    def generate_hub_file(self, shortLabel=None, longLabel=None, email=None):

        # set defaults unless given
        if shortLabel is None:
            shortLabel = self.genome

        if longLabel is None:
            longLabel = "Track hub for {}".format(self.genome)

        if email is None:
            email = "kgmcchesney@wisc.edu"

        # open the hub file
        hub_path = os.path.join(self.root_dir,"hub.txt")
        hub_file = open(hub_path,"w")

        # write it
        hub_file.write('''\
hub {}
shortLabel {}
longLabel {}
genomesFile genomes.txt
email {}
descriptionUrl about.html
'''.format(self.genome, shortLabel, longLabel, email))

        hub_file.close()

    def generate_genomes_file(self, length, description=None, organism=None, order=None, sci_name=None):
        
        if description is None:
            description = "{} Track Hub".format(self.genome)

        if organism is None:
            organism = self.genome

        if order is None:
            order = 1

        if sci_name is None:
            sci_name = self.genome

        # open the hub file
        genome_path = os.path.join(self.root_dir,"genomes.txt")
        genome_file = open(genome_path,"w")

        # write it
        genome_file.write('''\
genome {genome}
trackDb {genome}/trackDb.txt
groups {genome}/groups.txt
description {0}
twoBitPath {genome}/{genome}.2bit
organism {1}
defaultPos {genome}:0-{2}
orderKey {3}
scientificName {4}
htmlPath {genome}/description.html
'''.format(description, organism, length, order, sci_name, genome=self.genome))

        genome_file.close()

    # names and groups required!
    def generate_track_db_file(self, file_names, groups, long_labels=None):

        # create the file
        # needs to be in GENOME dir
        track_db = open(os.path.join(self.genome_dir,"trackDb.txt"), "a+")

        # groups and long labels and file_names must equal
        if not len(file_names) == len(groups):
            print "You must pass an equal number of groups and file names"
            raise SystemExit

        if long_labels is not None and len(long_labels) != len(file_names):
            print "You must pass an equal number of long labels and file names"
            raise SystemExit

        index = 0
        for file in file_names:
            name, ext = os.path.splitext(file)
            long_label = ""
            if long_labels is None:
                long_label = "{} track over {}".format(name, self.genome)

            else:
                long_label = long_labels[index]

            # big bed
            if ext == ".bb":
                track_db.write('''\
track {name}
longLabel {long_label}
shortLabel {name}
priority 1
visibility pack
group {group}
bigDataUrl {file}
type bigBed 9
itemRgb On\n
'''.format(name=name, long_label=long_label, group=groups[index], file=file))

            # big wig
            elif ext == ".bw":
                track_db.write('''\
track {name}
longLabel {long_label}
shortLabel {name}
priority 1
visibility full
group {group}
bigDataUrl {file}
type bigWig
maxHeightPixels 100:100:100\n
'''.format(name=name, long_label=long_label, group=groups[index], file=file))             

            # updated pointer
            index += 1

        # out of loop
        track_db.close()

    def generate_groups_file(self,names,labels):
        
        # check lengths
        if len(names) != len(labels):
            print "Groups names must be same length as group labels"
            raise SystemExit

        # open the file
        genomes_file = open(os.path.join(self.genome_dir,"genomes.txt"),"a+")

        priority = 1
        for name,label in zip(names,labels):
            genomes_file.write('''\
name {name}
label {label}
priority {priority}
defaultIsClosed 0\n
'''.format(name=name, label=label, priority=priority))
            priority += 1

        # close er up
        genomes_file.close()

