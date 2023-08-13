#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 06:40:21 2023

@author: hagar
"""

import pandas as pd
from pipelines.generalPipeline import general_pipe, DEPTH
import subprocess
from os import listdir
import glob
from utils.utils import get_sequences
import os
import pysam
from statistics import mean
import csv

SCRIPT_PATH = os.path.dirname(__file__)
cat = "cat %(cns)s* > %(aln_path)sall_not_aligned.fasta"
align = "augur align --sequences %(not_aligned)s --reference-sequence %(reference)s --output %(aligned)s"
#*gene_regions*#
pr_reg = (2252,2550)
rt_reg = (2661,3294)
int_reg = (4230, 5094)

class hiv(general_pipe):
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq)    
 
    
    def write_sub_fasta(self, fasta, path, regions, gene):
        with open(path + gene + ".fasta",'w') as f:
            start = regions[0]
            end = regions[1]
            for header, seq in fasta.items():
                f.write(">" + header + '\n')
                f.write(seq[start-1:end-1] + '\n')    
        
    def cut_to_gene(self, cns_path, aln_path):
        subprocess.call(cat % dict(cns=cns_path, aln_path=aln_path), shell=True)
        subprocess.call(align % dict(not_aligned=aln_path+"all_not_aligned.fasta", reference = self.reference, aligned=aln_path+"all_aligned.fasta"), shell=True)
        fasta = get_sequences(aln_path + "all_aligned.fasta")
        self.write_sub_fasta(fasta, cns_path, pr_reg, "reg_PR")
        self.write_sub_fasta(fasta, cns_path, rt_reg, "reg_RT")
        self.write_sub_fasta(fasta, cns_path, int_reg, "reg_Int")
        
        
        
    def excel_fasta(self, cns_path, forma):
        '''
        generate a report of each sample and gene(region_protein) fasta. merge it with MAGIC format. HIV department suppose to provid the MAGIC format

        Parameters
        ----------
        cns_path : consensus fasta path
        forma : MAGIC fomat excel file

        -------

        '''
        df = pd.DataFrame(columns=["SAMPLE_No_NGS", "Region_Protein", "fasta", "%coverage"])
        for file in listdir(cns_path):
            if "reg" in file:
                seqs = get_sequences(cns_path + file)
                gene = file.split("reg_")[1].split(".")[0]
                for sample, fasta in seqs.items():
                    if not sample == "K03455.1":
                        line = pd.DataFrame(data = {"SAMPLE_No_NGS" : sample,
                                                    "Region_Protein" : gene,
                                                    'fasta': fasta,
                                                    "%coverage": round(100 - fasta.count('N')/len(fasta))}
                                            , index=[0])
                        df = pd.concat([df,line], axis = 0)
        df.to_csv("QC/fasta.csv", index=False)
        forma = pd.read_excel(forma)
        mergi = pd.merge(forma, df, on=["SAMPLE_No_NGS","Region_Protein"], how = "left")
        mergi.to_excel("QC/final_report.xlsx", index=False)
        

    #@override
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path, cnsThresh):
        
    
        
        vcf_path = bam_path.replace("BAM","VCF")
        for file in os.listdir(vcf_path):
            if file.endswith("csv"):
                vcf = pd.read_csv(vcf_path + file)
                sample = file.split(".csv")[0]
                vcf["base"] = vcf[["%A","%T","%C","%G"]].apply(lambda row: row[row > 20].nlargest(2).index.values, axis=1)
                #clean
                vcf["base"] = vcf["base"].astype(str).apply(lambda row: row.replace("%", "").replace("'","").replace("[","").replace("]",""))
                vcf.loc[vcf["base"] == "", "base"] = 'A C T G'
                deg_nuc = pd.read_csv(SCRIPT_PATH + "/degenerate_nuc.csv", sep='\t')
                
                vcf = pd.merge(vcf, deg_nuc, on=["base"], how = "left")
                vcf.to_csv(vcf_path + sample + ".csv", index = False)
                
                cns = "".join(vcf["cns"].to_list())
               
                with open(cns_path +sample + ".fasta", 'w') as f:
                    f.write(">" + sample + '\n' + cns + '\n')
            
        for bam_file in os.listdir(bam_path): 
            if "sorted" in bam_file:
                sample = bam_file.split(".bam")[0].replace(".mapped.sorted", "")
                subprocess.call(DEPTH % dict(bam_path=bam_path, bam_file=bam_file, depth_path=depth_path, sample=sample), shell=True) 
            
    def results_report(self, bam_path, depth_path, output_report):
        
        fasta = pd.read_csv(output_report.replace("report", "fasta.csv"))
        
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
        writer.writerow(['sample', '%mapped','mapped_reads','total_reads','cov_bases','%coverage_RT','%coverage_PR',"%coverage_INT", 'mean_depth','max_depth','min_depth'])
        
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    total_reads = pysam.AlignmentFile(bam_path+bam_file.split(".mapped")[0]+".bam").count(until_eof=True) #need to fix 
                    cover_rt = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "RT"), "%coverage"].reset_index(drop=True)[0]
                    cover_pr = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "PR"), "%coverage"].reset_index(drop=True)[0]
                    cover_int = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "Int"), "%coverage"].reset_index(drop=True)[0]
                    
                    coverage_stats = pysam.coverage(bam_path+bam_file).split("\t")
                    mapped_reads = int(coverage_stats[11])
                    mapped_percentage = round(mapped_reads/total_reads*100,4) if total_reads else ''
                    cov_bases =  int(coverage_stats[12])
            
                    #depth 
                    depths = [int(x.split('\t')[2]) for x in open(depth_path+sample+".txt").readlines()]
                    depths = [i for i in depths if i != 0]
                    mean_depth = str(round(mean(depths),3)) if depths else ''
                    min_depth = min(depths) if depths else ''
                    max_depth = max(depths) if depths else ''
                    
                                        
                    writer.writerow([sample, mapped_percentage, mapped_reads, total_reads, cov_bases, cover_rt, cover_pr, cover_int, mean_depth, max_depth, min_depth])
                
        f.close()