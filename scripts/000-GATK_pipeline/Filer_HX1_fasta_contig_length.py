#coding:utf-8


import os
import sys

ref_fasta = sys.argv[1]
ref_fasta_fai = sys.argv[2]
fout_fasta = sys.argv[3]

ref_fasta = "/picb/humpopg-bigdata/xiebo/NGS/EAGLE_test/reference/HX1_fasta/hx1f4s4_3rdfixedv2.fa"
ref_fasta_fai = "/picb/humpopg-bigdata/xiebo/NGS/EAGLE_test/reference/HX1_fasta/hx1f4s4_3rdfixedv2.fa.fai"

ref_chr_list = []
with open(ref_fasta_fai) as file:
    for line in file:
        j = line.strip().split()
        chr = j[0]
        len = j[1]
        if int(len) > 1000:
            ref_chr_list.append(chr)

print(ref_chr_list[0:5])

fout_fasta = "/picb/humpopg-bigdata/xiebo/NGS/EAGLE_test/reference/HX1_fasta/hx1f4s4_3rdfixedv2.filter_len.fa"
if os.path.exists(fout_fasta):
    os.system("rm %s" %(fout_fasta))
else:
    print("start")
    os.system("samtools faidx %s %s > %s" %(ref_fasta,ref_chr_list[0],fout_fasta))
    for chr in ref_chr_list[1:]:
        os.system("samtools faidx %s %s >> %s" %(ref_fasta,chr,fout_fasta))
