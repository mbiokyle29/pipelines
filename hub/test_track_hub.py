from track_hub import TrackHub
import unittest
import os

class TestTrackHubInitalization(unittest.TestCase):
	dir_name = "./NEW-DIR/"
	genome = "EBV"
	shortLabel = "TEST_HUB"
	longLabel = "Long Label Test Hub"
	email = "test@email.com"
	length = 171823
	description = "TEST-DESC"
	organism = "EBV1"
	sci_name = "Human Herpes Virus 4"
	order = 1

	def tearDown(self):
		inner = os.path.join(self.dir_name, self.genome)
		hub_file = os.path.join(self.dir_name, "hub.txt")
		genome_file = os.path.join(self.dir_name, "genomes.txt")
		track_db_file = os.path.join(self.dir_name, self.genome, "trackDb.txt")

		for file in [hub_file, genome_file]:#, track_db_file]:
			if os.path.exists(file):
				os.unlink(file)
		
		os.rmdir(inner)
		os.rmdir(self.dir_name)

	def test_creation_with_new_dir(self):
		hub = TrackHub(self.dir_name, self.genome)

		self.assertEqual(hub.root_dir, self.dir_name)
		self.assertTrue(os.path.isdir(hub.root_dir))

	def test_creation_with_existing(self):
		os.mkdir(self.dir_name)
		hub = TrackHub(self.dir_name, self.genome)
		self.assertTrue(os.path.isdir(hub.root_dir))

	def test_generate_hub_file_with_params(self):
		hub = TrackHub(self.dir_name, self.genome)
		hub.generate_hub_file(shortLabel=self.shortLabel, longLabel=self.longLabel, email=self.email)
		hub_path = os.path.join(hub.root_dir,"hub.txt")

		# check if it was created
		self.assertTrue(os.path.exists(hub_path))

		# make sure contents match
		expected_body ='''\
hub {}
shortLabel {}
longLabel {}
genomesFile genomes.txt
email {}
descriptionUrl about.html
'''.format(self.genome, self.shortLabel, self.longLabel, self.email)
		hub_contents = open(hub_path).read()
		self.assertEqual(expected_body, hub_contents)

	def test_generate_hub_file_with_def(self):
		hub = TrackHub(self.dir_name, self.genome)
		hub.generate_hub_file()
		hub_path = os.path.join(hub.root_dir,"hub.txt")

		# check if it was created
		self.assertTrue(os.path.exists(hub_path))
		
		# make sure contents match
		expected_body ='''\
hub {}
shortLabel {}
longLabel Track hub for {}
genomesFile genomes.txt
email kgmcchesney@wisc.edu
descriptionUrl about.html
'''.format(self.genome, self.genome, self.genome, self.email)
		hub_contents = open(hub_path).read()
		self.assertEqual(expected_body, hub_contents)

	def test_generate_genome_file_with_params(self):
		hub = TrackHub(self.dir_name, self.genome)
		
		hub.generate_genomes_file(self.length, description=self.description, 
			organism=self.organism, order=self.order, sci_name=self.sci_name)

		genome_path = os.path.join(hub.root_dir,"genomes.txt")

		# check if it was created
		self.assertTrue(os.path.exists(genome_path))

		# make sure contents match
		expected_body ='''\
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
'''.format(self.description, self.organism, self.length, self.order, self.sci_name, genome=self.genome)
		genome_contents = open(genome_path).read()
		self.assertEqual(expected_body, genome_contents)

	def test_generate_track_db(self):
		hub = TrackHub(self.dir_name, self.genome)
		files = ["test.bb", "test.bw"]
		groups = ["geneAnnot", "RNAseq"]
		long_labels = ["EBV GENES BED FILE", "REP1 RNA from MeDIP"]
		hub.generate_track_db_file(files, groups, long_labels=long_labels)

		self.assertTrue(os.path.exists(hub.genome_dir+"/trackDb.txt"))

if __name__ == '__main__':
    unittest.main()