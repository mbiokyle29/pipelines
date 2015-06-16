rm */*.md5
rm */*.html
rm */*.htm
mv */*.fastq ./
rmdir AH5HWCBCXX AH5HYHBCXX BH5J2VBCXX
rename 's/_L00\d_R1_001//' *.gz
rename 's/^(\d)(.*)(\.fastq\.gz)/$2-rep$1$3/' *.gz
rename 's/IAhmC/CaFBS-AKATA-hMet/' *.gz
rename 's/IAmC/CaFBS-AKATA-Met/' *.gz
rename 's/UAhmC/UI-AKATA-hMet/' *.gz
rename 's/UAmC/UI-AKATA-Met/' *.gz
rename 's/IAI/CaFBS-AKATA-Input/' *.gz
rename 's/IhmC/CaFBS-hMet/' *.gz
rename 's/UAI/UI-AKATA-Input/' *.gz
rename 's/UhmC/UI-hMet/' *.gz
rename 's/II/CaFBS-Input/' *.gz
rename 's/ImC/CaFBS-Met/' *.gz
rename 's/UI/UI-Input/' *.gz
rename 's/UmC/UI-Met/' *.gz
