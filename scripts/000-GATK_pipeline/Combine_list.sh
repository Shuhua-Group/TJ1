pwd='/path/to/HX1/'
cd $pwd
mkdir list; cd list

#######prepare gvcf 
##Tujia HAN
ref_fai_path="/path/to/GWHAAAS00000000.genome.fasta.fai"
ref_chr_list=`cat $ref_fai_path | awk '{print $1}'`
for chr in $ref_chr_list
do
for i in HX1
do 
echo "path/to/${i}.${chr}.g.vcf.gz" >> ${chr}.list
done
done
