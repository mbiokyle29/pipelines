Eric's Workflow:

```
bwa index all-HPVs.fasta

cat run327.10_GAGTGG_L001_R1.fastq run327.10_GAGTGG_L002_R1.fastq run327.10_GAGTGG_L003_R1.fastq run327.10_GAGTGG_L004_R1.fastq > run327-10C.fastq

mkdir run327-10C.d;
cd run327-10C.d;
bwa aln ../all-HPVs.fasta ../run327-10C.fastq > aln_all-HPVs.sai;
bwa samse ../all-HPVs.fasta aln_all-HPVs.sai ../run327-10C.fastq > aln_all-HPVs.sam;
samtools view -F 4 -bT ../all-HPVs.fasta aln_all-HPVs.sam > aln_all-HPVs.bam;
samtools sort aln_all-HPVs.bam aln_all-HPVs-sorted;
samtools view aln_all-HPVs-sorted.bam > aln_all-HPVs-sorted.sam

***** stop here  ******

this work flow automatically extracts individual virus reads into separate files:

awk '$3 ~ /^EBV$/' aln_all-HPVs-sorted.sam > run327-10C_EBV.sam;
awk '$3 ~ /^HPV16$/' aln_all-HPVs-sorted.sam > run327-10C_HPV16.sam;
awk '$3 ~ /^HPV18$/' aln_all-HPVs-sorted.sam > run327-10C_HPV18.sam;
count run327-10C_EBV.sam; count run327-10C_HPV16.sam; count run327-10C_HPV18.sam; python ../Sam2Wig-BFR.py run327-10C_HPV16.sam ../HPV16-4k.fasta; python ../Sam2Wig-BFR.py run327-10C_EBV.sam ../EBV.fasta


rm aln_viruses.sai
rm aln_viruses.sam
rm aln_viruses.bam
cd ..
rm run327-10C.fastq
```

# NPC
SRR1654790
SRR1654791
SRR1654792
SRR1654793