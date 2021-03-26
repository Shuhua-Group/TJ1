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

## Introduction of pipeline
* Filter_data_and_use_rtg-tools.py: Filter variants and use rtg software to get the accuracy of variant calling  
* TJ1_SV.circlize.R: Use circlize packages to plot TJ1 SV distribution  
* run_EAGLE.sh: Use EAGLE to simulate a sample and its sequencing
