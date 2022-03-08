# TJ1

## Description
We applied multiple sequencing technologies to de novo assemble an individual genome (TJ1) from the Tujia population, an ethnic minority group most closely related to the Han Chinese. We used the following softwares to perform the variants calling, genotying and phasing with TJ1 as reference genome, and evaluate the performance with customized scripts.

## Softwares:
* EAGLE, https://github.com/sequencing/EAGLE  
* GATK, https://gatk.broadinstitute.org/hc/en-us  
* minimap2, https://github.com/lh3/minimap2  
* vg, https://github.com/vgteam/vg
* paragraph, https://github.com/Illumina/paragraph
* RTG tools, https://github.com/RealTimeGenomics/rtg-tools
* SHAPEIT4, https://odelaneau.github.io/shapeit4/  

## Introduction of scripts
* 000-GATK_pipeline: Short reads mapping and variants calling
* Filter_data_and_use_rtg-tools.py: Filter variants and use rtg software to get the accuracy of variant calling  
* TJ1_SV.circlize.R: Use circlize packages to plot TJ1 SV distribution  
* run_EAGLE.sh: Use EAGLE to simulate a sample and its sequencing

## Citation
Lou H, Gao Y, Xie B, Wang Y, Zhang H, Shi M, Ma S, Zhang X, Liu C, Xu S. Haplotype-resolved de novo assembly of a Tujia genome suggests the necessity for high-quality population-specific genome references. Under Review.
