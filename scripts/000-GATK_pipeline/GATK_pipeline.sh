#1. Input fasta file of the HX1 sample, filter out the contigs less than 1Kb in length, and then use the EAGLE software to simulate the short reads.
ref_HX1_fa=hx1f4s4_3rdfixedv2.fa
ref_HX1_fai=hx1f4s4_3rdfixedv2.fa.fai
ref_NH1_fa=GWHAAAS00000000.genome.fasta
ref_NH1_fai=GWHAAAS00000000.genome.fasta.fai
ref_HX1_filtered_fa=hx1f4s4_3rdfixedv2.filtered.fa

python Filter_HX1_fasta_contig_length.py ${ref_HX1_fa} ${ref_HX1_fai} ${ref_HX1_filtered_fa}

path/to/configureEAGLE.pl \
  --run-info=RunInfo_PairedReads2x251Cycles2x64Tiles.xml  \
  --reference-genome=${ref_HX1_filtered_fa} \
  --coverage-depth=30 \
  --motif-quality-drop-table=MotifQualityDropTables/DefaultMotifQualityDropTable.tsv \
  --template-length-table=TemplateLengthTables/TemplateLengthTableFrom2x250Run.tsv \
  --quality-table=NewQualityTable.read1.length251.qtable2 \
  --quality-table=NewQualityTable.read2.length251.qtable2
cd EAGLE
make fastq -j 15



#2. Use bwa software to map the simulated HX1 short reads to the NH1 reference genome, and then use the GATK process for variants calling.

#a. Mapping HX1 short reads to NH1 genome and remove duplicate reads
fastq_1=EAGLE_S1_L001_R1_001.fastq.gz
fastq_2=EAGLE_S1_L001_R2_001.fastq.gz

time bwa mem -M -t 10 -R "@RG\tID:HX1\tSM:HX1\tLB:HX1\tPU:HX1\tPL:ILLUMINA" \
   ${ref_NH1_fa} ${fastq_1} ${fastq_2}  | samtools view -bS - > HX1.pe.bam
time samtools sort -@ 5 -m 4G HX1.refNH1.pe.bam NH1; samtools index HX1.bam
time java -jar -Xmx4g -Djava.io.tmpdir=HX1/ MarkDuplicates.jar \
    INPUT=HX1.bam OUTPUT=HX1.dedup.bam \
    VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true METRICS_FILE=NH1.txt ASSUME_SORTED=true CREATE_INDEX=true
rm HX1.bam;rm HX1.pe.bam


#b. Using GATK software to do variant calling
ref_fai_path=${ref_NH1_fai}
ref_chr_list=`cat $ref_fai_path | awk '{print $1}'`
for chr in $ref_chr_list
    do
        time gatk --java-options \"-Xmx3G -XX:ParallelGCThreads=2 -Dsamjdk.compression_level=5 \" HaplotypeCaller \
        -R ${ref_NH1_fa} -ploidy 1 -L $chr -I HX1.dedup.bam -O HX1.$chr.g.vcf.gz -ERC GVCF -G StandardAnnotation \
        -G AS_StandardAnnotation -G StandardHCAnnotation --seconds-between-progress-updates 30
    done
sh Combine_list.sh;sh 170.JointCalling.sh
for chr in $ref_chr_list
    do
        sh 170.${chr}.sh
    done
sh 170.combine.sh


#c. GATK hard-filtering 
time gatk SelectVariants -select-type SNP -V HX1.genomewide.hc.vcf.gz -O HX1.genomewide.hc.snp.vcf.gz
time gatk VariantFiltration -V HX1.genomewide.hc.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" \
    --filter-name "Filter" -O HX1.genomewide.hc.snp.filter.vcf.gz
time gatk SelectVariants -select-type INDEL -V HX1.genomewide.hc.vcf.gz -O HX1.genomewide.hc.indel.vcf.gz
time gatk VariantFiltration -V HX1.genomewide.hc.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
    --filter-name "Filter" -O HX1.genomewide.hc.indel.filter.vcf.gz
time gatk MergeVcfs -I HX1.genomewide.hc.snp.filter.vcf.gz -I HX1.genomewide.hc.indel.filter.vcf.gz \
    -O HX1.genomewide.hc.filter.vcf.gz
vcftools --gzvcf HX1.genomewide.hc.filter.vcf.gz --remove-filtered-all --recode -c |bgzip -c > HX1.genomewide.hc.filtered.vcf.gz
