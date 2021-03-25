singularity exec -B /picb:/picb /picb/humpopg-bigdata5/xiebo/picb_eagle_2.5.1.sif configureEAGLE.pl \
  --run-info=/picb/humpopg-bigdata/xiebo/NGS/tools/EAGLE/data/RunInfo/RunInfo_PairedReads2x151Cycles1x1Tiles.xml  \
  --reference-genome=/picb/humpopg-bigdata5/xiebo/SV_calling/result/Pilon_test/All_polish_result/p0.NGS_350bp_8K.polish_3.fasta \
  --variant-list=/picb/humpopg-bigdata/xiebo/NGS/tools/EAGLE/data/Variants/None.vcf \
  --quality-table=/picb/humpopg-bigdata/xiebo/NGS/tools/EAGLE/data/QualityTables/tmp_1/qualityTable0.length151.read1.new.qval \
  --quality-table=/picb/humpopg-bigdata/xiebo/NGS/tools/EAGLE/data/QualityTables/tmp_1/qualityTable0.length151.read2.new.qval \
  --coverage-depth=30 
singularity exec -B /picb:/picb /picb/humpopg-bigdata5/xiebo/picb_eagle_2.5.1.sif make fastq -j 15
singularity exec -B /picb:/picb /picb/humpopg-bigdata5/xiebo/picb_eagle_2.5.1.sif make bam -j 15