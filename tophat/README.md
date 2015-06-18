
Here is erics workflow


cat run327.10_GAGTGG_L001_R1.fastq run327.10_GAGTGG_L002_R1.fastq run327.10_GAGTGG_L003_R1.fastq run327.10_GAGTGG_L004_R1.fastq > run327-10.fastq

bwa aln ../viruses.fasta run327-10.fastq > aln_viruses.sai; bwa samse ../viruses.fasta aln_viruses.sai run327-10.fastq > aln_viruses.sam

samtools view -F 4 -bT ../viruses.fasta aln_viruses.sam > aln_viruses.bam

samtools sort aln_viruses.bam aln_viruses-sorted; samtools view aln_viruses-sorted.bam > aln_viruses-sorted.sam

awk '$3 ~ /^EBV$/' aln_viruses-sorted.sam > run327-10_EBV.sam; awk '$3 ~ /^HPV16$/' aln_viruses-sorted.sam > run327-10_HPV16.sam; awk '$3 ~ /^HPV18$/' aln_viruses-sorted.sam > run327-10_HPV18.sam

count run327-10_EBV.sam; count run327-10_HPV16.sam; count run327-10_HPV18.sam


python ../Sam2Wig-BFR.py run327-10_HPV16.sam ../HPV16-4k.fasta; python ../Sam2Wig-BFR.py run327-10_EBV.sam ../EBV.fasta


rm aln_viruses.sai; rm aln_viruses.sam; rm aln_viruses.bam; rm run327-10.fastq

Here are the files as split

kyle@alpha-helix:/data/projects/HPV/fastqs > lss
total 150G
drwxr-xr-x 2 kyle 4.0K Jun 18 15:32 .
drwxr-xr-x 5 kyle 4.0K Jun 18 13:44 ..
-rw------- 1 kyle 4.4G Jan 11  2013 run327.9_CGTACG_L003_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.9_CGTACG_L004_R1.fastq
-rw------- 1 kyle 3.3G Jan 11  2013 run327.10_GAGTGG_L001_R1.fastq
-rw------- 1 kyle 3.5G Jan 11  2013 run327.10_GAGTGG_L002_R1.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_GAGTGG_L004_R1.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_GAGTGG_L003_R1.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_GTTTCG_L002_R1.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.9_CGTACG_L002_R1.fastq
-rw------- 1 kyle 4.3G Jan 11  2013 run327.9_CGTACG_L001_R1.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_GTTTCG_L004_R1.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.8_GTTTCG_L001_R1.fastq
-rw------- 1 kyle 4.1G Jan 11  2013 run327.8_GTTTCG_L003_R1.fastq
-rw------- 1 kyle 4.8G Jan 11  2013 run327.5_GTCCGC_L002_R1.fastq
-rw------- 1 kyle 3.1G Jan 11  2013 run327.7_GTGGCC_L004_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.7_GTGGCC_L002_R1.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.7_GTGGCC_L003_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_GTGAAA_L004_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_GTGAAA_L002_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_GTGAAA_L001_R1.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.6_GTGAAA_L003_R1.fastq
-rw------- 1 kyle 2.8G Jan 11  2013 run327.7_GTGGCC_L001_R1.fastq
-rw------- 1 kyle 4.7G Jan 11  2013 run327.5_GTCCGC_L004_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_GTCCGC_L003_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.4_CCGTCC_L002_R1.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.4_CCGTCC_L004_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_GTCCGC_L001_R1.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_CCGTCC_L003_R1.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_CCGTCC_L001_R1.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.3_ATGTCA_L004_R1.fastq
-rw------- 1 kyle 3.8G Jan 10  2013 run327.2_AGTTCC_L003_R1.fastq
-rw------- 1 kyle 4.1G Jan 10  2013 run327.3_ATGTCA_L002_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_AGTTCC_L004_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_AGTTCC_L002_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_ATGTCA_L001_R1.fastq
-rw------- 1 kyle 3.7G Jan 10  2013 run327.2_AGTTCC_L001_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_ATGTCA_L003_R1.fastq
-rw------- 1 kyle 3.2G Jan 10  2013 run327.1_AGTCAA_L002_R1.fastq
-rw------- 1 kyle 3.0G Jan 10  2013 run327.1_AGTCAA_L003_R1.fastq
-rw------- 1 kyle 3.1G Jan 10  2013 run327.1_AGTCAA_L004_R1.fastq
-rw------- 1 kyle 2.9G Jan 10  2013 run327.1_AGTCAA_L001_R1.fastq



kyle@alpha-helix:/data/projects/HPV/fastqs > rename "s/_[ATCG]//" *.fastq
kyle@alpha-helix:/data/projects/HPV/fastqs > lss
total 150G
drwxr-xr-x 2 kyle 4.0K Jun 18 15:39 .
drwxr-xr-x 5 kyle 4.0K Jun 18 13:44 ..
-rw------- 1 kyle 4.4G Jan 11  2013 run327.9GTACG_L003_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.9GTACG_L004_R1.fastq
-rw------- 1 kyle 3.3G Jan 11  2013 run327.10AGTGG_L001_R1.fastq
-rw------- 1 kyle 3.5G Jan 11  2013 run327.10AGTGG_L002_R1.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10AGTGG_L004_R1.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10AGTGG_L003_R1.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8TTTCG_L002_R1.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.9GTACG_L002_R1.fastq
-rw------- 1 kyle 4.3G Jan 11  2013 run327.9GTACG_L001_R1.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8TTTCG_L004_R1.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.8TTTCG_L001_R1.fastq
-rw------- 1 kyle 4.1G Jan 11  2013 run327.8TTTCG_L003_R1.fastq
-rw------- 1 kyle 4.8G Jan 11  2013 run327.5TCCGC_L002_R1.fastq
-rw------- 1 kyle 3.1G Jan 11  2013 run327.7TGGCC_L004_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.7TGGCC_L002_R1.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.7TGGCC_L003_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6TGAAA_L004_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6TGAAA_L002_R1.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6TGAAA_L001_R1.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.6TGAAA_L003_R1.fastq
-rw------- 1 kyle 2.8G Jan 11  2013 run327.7TGGCC_L001_R1.fastq
-rw------- 1 kyle 4.7G Jan 11  2013 run327.5TCCGC_L004_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5TCCGC_L003_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.4CGTCC_L002_R1.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.4CGTCC_L004_R1.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5TCCGC_L001_R1.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4CGTCC_L003_R1.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4CGTCC_L001_R1.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.3TGTCA_L004_R1.fastq
-rw------- 1 kyle 3.8G Jan 10  2013 run327.2GTTCC_L003_R1.fastq
-rw------- 1 kyle 4.1G Jan 10  2013 run327.3TGTCA_L002_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2GTTCC_L004_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2GTTCC_L002_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3TGTCA_L001_R1.fastq
-rw------- 1 kyle 3.7G Jan 10  2013 run327.2GTTCC_L001_R1.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3TGTCA_L003_R1.fastq
-rw------- 1 kyle 3.2G Jan 10  2013 run327.1GTCAA_L002_R1.fastq
-rw------- 1 kyle 3.0G Jan 10  2013 run327.1GTCAA_L003_R1.fastq
-rw------- 1 kyle 3.1G Jan 10  2013 run327.1GTCAA_L004_R1.fastq
-rw------- 1 kyle 2.9G Jan 10  2013 run327.1GTCAA_L001_R1.fastq


kyle@alpha-helix:/data/projects/HPV/fastqs > rename -n "s/_R1//" *.fastq
rename(run327.10AGTGG_L001_R1.fastq, run327.10AGTGG_L001.fastq)
rename(run327.10AGTGG_L002_R1.fastq, run327.10AGTGG_L002.fastq)
rename(run327.10AGTGG_L003_R1.fastq, run327.10AGTGG_L003.fastq)
rename(run327.10AGTGG_L004_R1.fastq, run327.10AGTGG_L004.fastq)
rename(run327.1GTCAA_L001_R1.fastq, run327.1GTCAA_L001.fastq)
rename(run327.1GTCAA_L002_R1.fastq, run327.1GTCAA_L002.fastq)
rename(run327.1GTCAA_L003_R1.fastq, run327.1GTCAA_L003.fastq)
rename(run327.1GTCAA_L004_R1.fastq, run327.1GTCAA_L004.fastq)
rename(run327.2GTTCC_L001_R1.fastq, run327.2GTTCC_L001.fastq)
rename(run327.2GTTCC_L002_R1.fastq, run327.2GTTCC_L002.fastq)
rename(run327.2GTTCC_L003_R1.fastq, run327.2GTTCC_L003.fastq)
rename(run327.2GTTCC_L004_R1.fastq, run327.2GTTCC_L004.fastq)
rename(run327.3TGTCA_L001_R1.fastq, run327.3TGTCA_L001.fastq)
rename(run327.3TGTCA_L002_R1.fastq, run327.3TGTCA_L002.fastq)
rename(run327.3TGTCA_L003_R1.fastq, run327.3TGTCA_L003.fastq)
rename(run327.3TGTCA_L004_R1.fastq, run327.3TGTCA_L004.fastq)
rename(run327.4CGTCC_L001_R1.fastq, run327.4CGTCC_L001.fastq)
rename(run327.4CGTCC_L002_R1.fastq, run327.4CGTCC_L002.fastq)
rename(run327.4CGTCC_L003_R1.fastq, run327.4CGTCC_L003.fastq)
rename(run327.4CGTCC_L004_R1.fastq, run327.4CGTCC_L004.fastq)
rename(run327.5TCCGC_L001_R1.fastq, run327.5TCCGC_L001.fastq)
rename(run327.5TCCGC_L002_R1.fastq, run327.5TCCGC_L002.fastq)
rename(run327.5TCCGC_L003_R1.fastq, run327.5TCCGC_L003.fastq)
rename(run327.5TCCGC_L004_R1.fastq, run327.5TCCGC_L004.fastq)
rename(run327.6TGAAA_L001_R1.fastq, run327.6TGAAA_L001.fastq)
rename(run327.6TGAAA_L002_R1.fastq, run327.6TGAAA_L002.fastq)
rename(run327.6TGAAA_L003_R1.fastq, run327.6TGAAA_L003.fastq)
rename(run327.6TGAAA_L004_R1.fastq, run327.6TGAAA_L004.fastq)
rename(run327.7TGGCC_L001_R1.fastq, run327.7TGGCC_L001.fastq)
rename(run327.7TGGCC_L002_R1.fastq, run327.7TGGCC_L002.fastq)
rename(run327.7TGGCC_L003_R1.fastq, run327.7TGGCC_L003.fastq)
rename(run327.7TGGCC_L004_R1.fastq, run327.7TGGCC_L004.fastq)
rename(run327.8TTTCG_L001_R1.fastq, run327.8TTTCG_L001.fastq)
rename(run327.8TTTCG_L002_R1.fastq, run327.8TTTCG_L002.fastq)
rename(run327.8TTTCG_L003_R1.fastq, run327.8TTTCG_L003.fastq)
rename(run327.8TTTCG_L004_R1.fastq, run327.8TTTCG_L004.fastq)
rename(run327.9GTACG_L001_R1.fastq, run327.9GTACG_L001.fastq)
rename(run327.9GTACG_L002_R1.fastq, run327.9GTACG_L002.fastq)
rename(run327.9GTACG_L003_R1.fastq, run327.9GTACG_L003.fastq)
rename(run327.9GTACG_L004_R1.fastq, run327.9GTACG_L004.fastq)

kyle@alpha-helix:/data/projects/HPV/fastqs > rename "s/[ATCG]+//" *.fastq
kyle@alpha-helix:/data/projects/HPV/fastqs > ls
run327.10_L001.fastq  run327.2_L001.fastq  run327.4_L001.fastq  run327.6_L001.fastq  run327.8_L001.fastq
run327.10_L002.fastq  run327.2_L002.fastq  run327.4_L002.fastq  run327.6_L002.fastq  run327.8_L002.fastq
run327.10_L003.fastq  run327.2_L003.fastq  run327.4_L003.fastq  run327.6_L003.fastq  run327.8_L003.fastq
run327.10_L004.fastq  run327.2_L004.fastq  run327.4_L004.fastq  run327.6_L004.fastq  run327.8_L004.fastq
run327.1_L001.fastq   run327.3_L001.fastq  run327.5_L001.fastq  run327.7_L001.fastq  run327.9_L001.fastq
run327.1_L002.fastq   run327.3_L002.fastq  run327.5_L002.fastq  run327.7_L002.fastq  run327.9_L002.fastq
run327.1_L003.fastq   run327.3_L003.fastq  run327.5_L003.fastq  run327.7_L003.fastq  run327.9_L003.fastq
run327.1_L004.fastq   run327.3_L004.fastq  run327.5_L004.fastq  run327.7_L004.fastq  run327.9_L004.fastq
kyle@alpha-helix:/data/projects/HPV/fastqs > lss
total 150G
drwxr-xr-x 2 kyle 4.0K Jun 18 15:40 .
drwxr-xr-x 5 kyle 4.0K Jun 18 13:44 ..
-rw------- 1 kyle 4.4G Jan 11  2013 run327.9_L003.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.9_L004.fastq
-rw------- 1 kyle 3.3G Jan 11  2013 run327.10_L001.fastq
-rw------- 1 kyle 3.5G Jan 11  2013 run327.10_L002.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_L004.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_L003.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_L002.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.9_L002.fastq
-rw------- 1 kyle 4.3G Jan 11  2013 run327.9_L001.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_L004.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.8_L001.fastq
-rw------- 1 kyle 4.1G Jan 11  2013 run327.8_L003.fastq
-rw------- 1 kyle 4.8G Jan 11  2013 run327.5_L002.fastq
-rw------- 1 kyle 3.1G Jan 11  2013 run327.7_L004.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.7_L002.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.7_L003.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L004.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L002.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L001.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.6_L003.fastq
-rw------- 1 kyle 2.8G Jan 11  2013 run327.7_L001.fastq
-rw------- 1 kyle 4.7G Jan 11  2013 run327.5_L004.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_L003.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.4_L002.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.4_L004.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_L001.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_L003.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_L001.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.3_L004.fastq
-rw------- 1 kyle 3.8G Jan 10  2013 run327.2_L003.fastq
-rw------- 1 kyle 4.1G Jan 10  2013 run327.3_L002.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_L004.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_L002.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_L001.fastq
-rw------- 1 kyle 3.7G Jan 10  2013 run327.2_L001.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_L003.fastq
-rw------- 1 kyle 3.2G Jan 10  2013 run327.1_L002.fastq
-rw------- 1 kyle 3.0G Jan 10  2013 run327.1_L003.fastq
-rw------- 1 kyle 3.1G Jan 10  2013 run327.1_L004.fastq
-rw------- 1 kyle 2.9G Jan 10  2013 run327.1_L001.fastq


A smarter way to catem
kyle@alpha-helix:/data/projects/HPV/fastqs > for i in `seq 1 10`; do lss run327.$i\_L00*.fastq; echo "LINE";  done
-rw------- 1 kyle 3.2G Jan 10  2013 run327.1_L002.fastq
-rw------- 1 kyle 3.0G Jan 10  2013 run327.1_L003.fastq
-rw------- 1 kyle 3.1G Jan 10  2013 run327.1_L004.fastq
-rw------- 1 kyle 2.9G Jan 10  2013 run327.1_L001.fastq
LINE
-rw------- 1 kyle 3.8G Jan 10  2013 run327.2_L003.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_L004.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.2_L002.fastq
-rw------- 1 kyle 3.7G Jan 10  2013 run327.2_L001.fastq
LINE
-rw------- 1 kyle 4.0G Jan 11  2013 run327.3_L004.fastq
-rw------- 1 kyle 4.1G Jan 10  2013 run327.3_L002.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_L001.fastq
-rw------- 1 kyle 3.9G Jan 10  2013 run327.3_L003.fastq
LINE
-rw------- 1 kyle 4.5G Jan 11  2013 run327.4_L002.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.4_L004.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_L003.fastq
-rw------- 1 kyle 4.4G Jan 11  2013 run327.4_L001.fastq
LINE
-rw------- 1 kyle 4.8G Jan 11  2013 run327.5_L002.fastq
-rw------- 1 kyle 4.7G Jan 11  2013 run327.5_L004.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_L003.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.5_L001.fastq
LINE
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L004.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L002.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.6_L001.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.6_L003.fastq
LINE
-rw------- 1 kyle 3.1G Jan 11  2013 run327.7_L004.fastq
-rw------- 1 kyle 3.0G Jan 11  2013 run327.7_L002.fastq
-rw------- 1 kyle 2.9G Jan 11  2013 run327.7_L003.fastq
-rw------- 1 kyle 2.8G Jan 11  2013 run327.7_L001.fastq
LINE
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_L002.fastq
-rw------- 1 kyle 4.2G Jan 11  2013 run327.8_L004.fastq
-rw------- 1 kyle 4.0G Jan 11  2013 run327.8_L001.fastq
-rw------- 1 kyle 4.1G Jan 11  2013 run327.8_L003.fastq
LINE
-rw------- 1 kyle 4.4G Jan 11  2013 run327.9_L003.fastq
-rw------- 1 kyle 4.5G Jan 11  2013 run327.9_L004.fastq
-rw------- 1 kyle 4.6G Jan 11  2013 run327.9_L002.fastq
-rw------- 1 kyle 4.3G Jan 11  2013 run327.9_L001.fastq
LINE
-rw------- 1 kyle 3.3G Jan 11  2013 run327.10_L001.fastq
-rw------- 1 kyle 3.5G Jan 11  2013 run327.10_L002.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_L004.fastq
-rw------- 1 kyle 3.4G Jan 11  2013 run327.10_L003.fastq
LINE
kyle@al


kyle@alpha-helix:/data/projects/HPV/fastqs > for i in `seq 1 10`
> do
> cat run327.$i\_L00*.fastq >> run327.$i.fastq
> done
