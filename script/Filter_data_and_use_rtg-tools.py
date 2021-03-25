#coding: utf-8

import os
import sys
import gzip 


workdir = os.getcwd()

rtg = workdir + "/bin/rtg"
ref_sdf = workdir + "/rtg_test/reference/TJ1p1.sdf"
ref_fasta_fai = workdir + "/p1.NGS_350bp_8K.polish_3.fasta.fai"
ref_mask_path = workdir + "../Minimap2_result/"
ref_gap_path =  workdir + "/Gap_area_filter/"
ref_True_variant_path = workdir + "/TJ1_Ture_variants.vcf"
ref_EAGLE_snp_path = workdir + "/TJ1_p1.snp.vcf.gz"
ref_EAGLE_indel_path = workdir + "/TJ1_p1.indel.vcf.gz"


fout_base_snp_path = workdir + "/base_snp.vcf"
fout_query_snp_path = workdir + "/query_snp.vcf"
fout_base_snp_file = open(fout_base_snp_path,"w")
fout_query_snp_file = open(fout_query_snp_path,"w")

ref_scaffold_dt = {}
with open(ref_True_variant_path) as file:
    for line in file:
        j = line.strip().split()
        scaffold = j[0]
        if scaffold not in ref_scaffold_dt:
            ref_scaffold_dt[scaffold] = []

def WHOLE(ref_WHOLE_path):
    lnum = 0
    allgenes=[]
    for line in ref_WHOLE_path:
        val = line.strip()
        j = val.split()
        chr = j[0]
        start = int(j[1])
        end = int(j[2])
        WHOLE_list = [chr,start,end]
        allgenes.append(WHOLE_list)
    return(allgenes)
    
def ingene(pos,allgenes):
    pos = int(pos)
    Gene_list = []
    for i in range(len(allgenes)):
        ele = allgenes[i]
        if pos<=ele[2] and pos>=ele[1]:
            gene_list = [i]
            Gene_list.append(gene_list)
        if pos<ele[1]:
            break            
    return(Gene_list)

def Ingap(start,end,allgenes):
    start = int(start)
    end = int(end)
    Gene_list = []
    for i in range(len(allgenes)):
        ele = allgenes[i]
        if ele[1]<=start<=ele[2]:
            gene_list = [i]
            Gene_list.append(gene_list)
        elif start<=ele[1]<=end:
            gene_list = [i]
            Gene_list.append(gene_list)
        elif end < ele[1]:
            break   
    return(Gene_list)


####################################################################
##### SNP
####################################################################
####Ture snp
ref_True_variant_list = []
ref_True_variant_dt = {}
ref_chr_dt = {}
with open(ref_True_variant_path) as file:
    for line in file:
        val = line.strip()
        if val.count("#") != 0:
            fout_base_snp_file.write(line)
        else:
            j = val.split("\t")
            chr = j[0]
            pos = j[1]
            start = str(int(pos) - 1)
            end = pos
            ref = j[3]
            alt = j[4]
            cbind = chr + '_' + start + "_" + end
            if (len(ref) == 1) and (len(alt) == 1):
                cbind = chr + "-" + pos + "-" + ref + "-" + alt
                if os.path.exists(ref_mask_path + chr + ".variant_cluster_area.txt"):
                    if chr not in ref_chr_dt:
                        ref_chr_dt[chr] = 0
                        ref_WHOLE_path = open(ref_mask_path + chr + ".variant_cluster_area.txt","r")
                        allgenes = WHOLE(ref_WHOLE_path)
                    Gene_list = ingene(pos,allgenes)
                    if len(Gene_list) != 0:
                        nu = Gene_list[0][0]
                        del allgenes[0:nu]
                    else:
                        fout_base_snp_file.write(line)
                        ref_True_variant_list.append(cbind)
                        ref_True_variant_dt[cbind] = 0
                else:
                    fout_base_snp_file.write(line)
                    ref_True_variant_list.append(cbind)
                    ref_True_variant_dt[cbind] = 0


#####EAGLE SNP
ref_EAGLE_variant_in_gap_list = []
ref_EAGLE_variant_list = []
ref_EAGLE_variant_dt = {}
ref_chr_dt = {}
ref_gap_chr_dt = {}
with gzip.open(ref_EAGLE_snp_path) as file:
    for line in file:
        val = line.strip()
        if val.count("#") != 0:
            fout_query_snp_file.write(line)
        else:
            j = val.split("\t")
            chr = j[0]
            if chr in ref_scaffold_dt:
                pos = j[1]
                start = str(int(pos) - 1)
                end = pos
                cbind = chr + '_' + start + "_" + end
                cluster_filter_nu = 0
                gap_filter_nu = 0
                ref = j[3]
                alt = j[4]
                FILTER = j[6]
                genotype = j[-1]
                cbind = chr + "-" + pos + "-" + ref + "-" + alt
                if os.path.exists(ref_mask_path + chr + ".variant_cluster_area.txt"):
                    if chr not in ref_chr_dt:
                        ref_chr_dt[chr] = 0
                        ref_WHOLE_path = open(ref_mask_path + chr + ".variant_cluster_area.txt","r")
                        allgenes = WHOLE(ref_WHOLE_path)
                    Gene_list = ingene(pos,allgenes)
                    if len(Gene_list) != 0:
                        nu = Gene_list[0][0]
                        del allgenes[0:nu]
                    else:
                        cluster_filter_nu = 1
                else:
                    cluster_filter_nu = 1

                if chr not in ref_gap_chr_dt:
                    ref_gap_chr_dt[chr] = 0
                    ref_WHOLE_path = open(ref_gap_path + chr + ".gap_area.txt")
                    allgaps = WHOLE(ref_WHOLE_path)
                    ref_WHOLE_path.close()
                Gap_list = Ingap(pos,pos,allgaps)
                if len(Gap_list) != 0:
                    nu = Gap_list[0][0]
                    del Gap_list[0:nu]
                else:
                    gap_filter_nu = 1
                if (cluster_filter_nu == 1) and (gap_filter_nu == 1):
                    fout_query_snp_file.write(line)
                    ref_EAGLE_variant_list.append(cbind)
                    ref_EAGLE_variant_dt[cbind] = 0
                if gap_filter_nu == 0:
                    ref_EAGLE_variant_in_gap_list.append(cbind)

fout_base_snp_file.close()
fout_query_snp_file.close()


os.system("bgzip base_snp.vcf")
os.system("tabix -p vcf base_snp.vcf.gz")    
os.system("bgzip query_snp.vcf")
os.system("tabix -p vcf query_snp.vcf.gz")    
os.system("%s vcfeval -b base_snp.vcf.gz -c query_snp.vcf.gz -o snp_output -t %s >/dev/null 2>&1" %(rtg,ref_sdf))


def count_snp_nu(ref_path):
    snp_nu = 0
    with gzip.open(ref_path) as file:
        for line in file:
            val = line.strip()
            if val.count("#") != 0:
                pass
            else:
                j = val.split("\t")
                chr = j[0]
                pos = j[1]
                ref = j[3]
                alt = j[4]
                if (len(ref) == 1) and (len(alt) == 1):
                    snp_nu += 1
    return(snp_nu)

ref_True_path = "./base_snp.vcf.gz"
ref_GATK_path = "./query_snp.vcf.gz"
ref_GATK_TP_path = "./snp_output/tp.vcf.gz"

True_snp = count_snp_nu(ref_True_path)
GATK_snp = count_snp_nu(ref_GATK_path)
GATK_TP_snp = count_snp_nu(ref_GATK_TP_path)

print('#############SNP##############')
print("Precision_rate=" + str(float(GATK_TP_snp)/int(GATK_snp)))
print("recall_rate=" + str(float(GATK_TP_snp)/int(True_snp)))
print("GATK_TP_snp_nu=" + str(GATK_TP_snp))
print("GATK_snp_nu=" + str(GATK_snp)) 
print("True_snp_nu=" + str(True_snp)) 


####################################################################
##### INDEL
####################################################################
######Ture INDEL
fout_base_indel_path = "base_indel.vcf"
fout_query_indel_path = "query_indel.vcf"
fout_base_indel_file = open(fout_base_indel_path,"w")
fout_query_indel_file = open(fout_query_indel_path,"w")


ref_True_variant_list = []
ref_True_variant_dt = {}
ref_chr_dt = {}
with open(ref_True_variant_path) as file:
    for line in file:
        val = line.strip()
        if val.count("#") != 0:
            fout_base_indel_file.write(line)
        else:
            j = val.split("\t")
            chr = j[0]
            pos = j[1]
            start = str(int(pos) - 1)
            end = pos
            ref = j[3]
            alt = j[4]
            cbind = chr + '_' + start + "_" + end
            if (len(ref) != 1) or (len(alt) != 1):
                if (len(ref) <= 50) and (len(alt) <= 50):
                    cbind = chr + "-" + pos + "-" + ref + "-" + alt
                    if os.path.exists(ref_mask_path + chr + ".variant_cluster_area.txt"):
                        if chr not in ref_chr_dt: 
                            ref_chr_dt[chr] = 0
                            ref_WHOLE_path = open(ref_mask_path + chr + ".variant_cluster_area.txt","r")
                            allgenes = WHOLE(ref_WHOLE_path)
                        Gene_list = ingene(pos,allgenes)
                        if len(Gene_list) != 0:
                            nu = Gene_list[0][0]
                            del allgenes[0:nu]
                        else:
                            fout_base_indel_file.write(line)
                            ref_True_variant_list.append(cbind)
                            ref_True_variant_dt[cbind] = 0
                    else:
                        fout_base_indel_file.write(line)
                        ref_True_variant_list.append(cbind)
                        ref_True_variant_dt[cbind] = 0

########EAGLE INDEL
ref_EAGLE_variant_list = []
ref_EAGLE_variant_dt = {}
ref_EAGLE_variant_in_gap_list = []
ref_chr_dt = {}
ref_gap_chr_dt = {}
with gzip.open(ref_EAGLE_indel_path) as file:
    for line in file:
        val = line.strip()
        if val.count("#") != 0:
            fout_query_indel_file.write(line)
        else:
            j = val.split("\t")
            chr = j[0]
            if chr in ref_scaffold_dt:
                cluster_filter_nu = 0
                gap_filter_nu = 0
                pos = j[1]
                ref = j[3]
                alt = j[4]
                start = str(int(pos) - 1)
                end = pos
                cbind = chr + '_' + start + "_" + end
                if (len(ref) <= 50) and (len(alt) <= 50):
                    FILTER = j[6]
                    genotype = j[-1]
                    cbind = chr + "-" + pos + "-" + ref + "-" + alt
                    if os.path.exists(ref_mask_path + chr + ".variant_cluster_area.txt"):
                        if chr not in ref_chr_dt:
                            ref_chr_dt[chr] = 0
                            ref_WHOLE_path = open(ref_mask_path + chr + ".variant_cluster_area.txt","r")
                            allgenes = WHOLE(ref_WHOLE_path)
                            ref_WHOLE_path.close()
                        Gene_list = ingene(pos,allgenes)
                        if len(Gene_list) != 0:
                            nu = Gene_list[0][0]
                            del allgenes[0:nu]
                        else:
                            cluster_filter_nu = 1
                    else:
                        cluster_filter_nu = 1
                    if chr not in ref_gap_chr_dt:
                        ref_gap_chr_dt[chr] = 0
                        ref_WHOLE_path = open(ref_gap_path + chr + ".gap_area.txt")
                        allgaps = WHOLE(ref_WHOLE_path)
                        ref_WHOLE_path.close()
                    end = str(int(pos) + abs(len(ref)-len(alt)))
                    Gap_list = Ingap(pos,end,allgaps)
                    if len(Gap_list) != 0:
                        nu = Gap_list[0][0]
                        del Gap_list[0:nu]
                    else:
                        gap_filter_nu = 1
                    if (cluster_filter_nu == 1) and (gap_filter_nu == 1):
                        fout_query_indel_file.write(line)
                        ref_EAGLE_variant_list.append(cbind)
                        ref_EAGLE_variant_dt[cbind] = 0
                    if gap_filter_nu == 0:
                        ref_EAGLE_variant_in_gap_list.append(line)

fout_base_indel_file.close()
fout_query_indel_file.close()
os.system("bgzip base_indel.vcf")
os.system("tabix -p vcf base_indel.vcf.gz")    
os.system("bgzip query_indel.vcf")
os.system("tabix -p vcf query_indel.vcf.gz")    
os.system("%s vcfeval -b base_indel.vcf.gz -c query_indel.vcf.gz -o indel_output -t %s >/dev/null 2>&1" %(rtg,ref_sdf))


def count_indel_nu(ref_path):
    indel_nu = 0
    with gzip.open(ref_path) as file:
        for line in file:
            val = line.strip()
            if val.count("#") != 0:
                pass
            else:
                j = val.split("\t")
                chr = j[0]
                pos = j[1]
                ref = j[3]
                alt = j[4]
                if (len(ref) != 1) or (len(alt) != 1):
                    if (len(ref) <= 50) and (len(alt) <= 50):
                        indel_nu += 1
    return(indel_nu)

ref_True_path = "./base_indel.vcf.gz"
ref_GATK_path = "./query_indel.vcf.gz"
ref_GATK_TP_path = "./indel_output/tp.vcf.gz"

True_indel = count_indel_nu(ref_True_path)
GATK_indel = count_indel_nu(ref_GATK_path)
GATK_TP_indel = count_indel_nu(ref_GATK_TP_path)


print('#############INDEL##############')
print("Precision_rate=" + str(float(GATK_TP_indel)/int(GATK_indel)))
print("recall_rate=" + str(float(GATK_TP_indel)/int(True_indel)))
print("GATK_TP_indel_nu=" + str(GATK_TP_indel))
print("GATK_indel_nu=" + str(GATK_indel)) 
print("True_indel_nu=" + str(True_indel))
