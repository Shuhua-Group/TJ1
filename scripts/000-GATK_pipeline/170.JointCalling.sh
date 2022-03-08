#!/bin/bash

cwd='path/to/HX1_refNH1/'

mkdir -p $cwd

###ref files and bins
fasta='path/to/GWHAAAS00000000.genome.fasta'

###############################################################
proid="NH1"
mkdir -p $cwd/170.GenotypeGVCFs.joint.calling; cd $cwd/170.GenotypeGVCFs.joint.calling
ref_fai_path="/path/to/GWHAAAS00000000.genome.fasta.fai"
ref_chr_list=`cat $ref_fai_path | awk '{print $1}'`
for chr in $ref_chr_list
do
	bash=170.$chr.sh
	echo "#$ -N genovcf.$chr" > $bash
	echo "mkdir -p tmp/${chr}" >> $bash 
	echo "time java -jar -Xmx20g -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${chr} path/to/GenomeAnalysisTK3.6.jar -T GenotypeGVCFs -R $fasta \\" >> $bash
	for line in `cat $cwd/list/$chr.list`
	do
		echo "--variant $line \\" >> $bash
	done
	echo "-L $chr -nt 8  -o $proid.$chr.hc.vcf.gz" >> $bash
done
#combine as one file
vlist=" "
ref_fai_path="path/to/GWHAAAS00000000.genome.fasta.fai"
ref_chr_list=`cat $ref_fai_path | awk '{print $1}'`
for chr in $ref_chr_list
do
	vlist=" $vlist -V "$proid.$chr.hc.vcf.gz
done
bash=170.combine.sh
echo "#$ -N com170" >$bash
echo "mkdir -p tmp/combine" >> $bash 
echo "time java -Xmx256g -Djava.io.tmpdir=tmp/combine -cp path/to/GenomeAnalysisTK3.6.jar  org.broadinstitute.gatk.tools.CatVariants  -R $fasta  $vlist  -out $proid.genomewide.hc.vcf.gz -assumeSorted"  >>$bash
