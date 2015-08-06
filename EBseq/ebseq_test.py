import unittest
import logging
import re
import ebseq_extras

class TestEbseqExtras(unittest.TestCase):

	def setUp(self):
		self.conf_file = "../ebseq-test/conf.tsv"
		self.log = logging
		self.extras = ebseq_extras.EbseqExtras(self.log)

	def testReadConf(self):

		# will throw error in test if fails
		self.extras.read_configuration(self.conf_file)

		# inspect the records to see if made right
		num_recs = 9
		self.assertEqual(len(self.extras.samples), num_recs)
		for rec in self.extras.samples:
			self.assertTrue(rec.mean in [260,100])
			self.assertTrue(re.match("[ABC]-rep[123]\.fastq",rec.name) is not None)

	def testGenSampleList(self):

		# read it in first
		self.extras.read_configuration(self.conf_file)

		sample_string = self.extras.gen_sample_list()
		expected  = "A-rep1.genes.results A-rep2.genes.results A-rep3.genes.results "
		expected += "B-rep1.genes.results B-rep2.genes.results B-rep3.genes.results "
		expected += "C-rep1.genes.results C-rep2.genes.results C-rep3.genes.results"

		self.assertEqual(sample_string, expected)

	def testMeanLookup(self):

		# read it in first
		self.extras.read_configuration(self.conf_file)
		test = "A-rep1.fastq"
		
		self.assertEqual(self.extras.get_mean_length(test), 100)

	def testCondStr(self):

		# read it in first
		self.extras.read_configuration(self.conf_file)
		expected = "3,3,3"
		
		self.assertEqual(self.extras.gen_cond_string(), expected)


if __name__ == '__main__':
    unittest.main()