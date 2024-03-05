#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 06:40:21 2023

@author: hagar
"""

import pandas as pd
from pipelines.generalPipeline import general_pipe, DEPTH, ALL_NOT_ALIGNED
import subprocess
from os import listdir
from utils.utils import get_sequences, write_sub_fasta, mafft
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
    def __init__(self, reference, fastq, threads, minion, metadata):
        super().__init__(reference, fastq, minion, threads)    
        self.metadata = metadata
 
    
    def cut_genes(self, aln_path):
        fasta = get_sequences(aln_path + "all_aligned.fasta")
        write_sub_fasta(fasta, aln_path, pr_reg, "reg_PR")
        write_sub_fasta(fasta, aln_path, int_reg, "reg_Int")
        write_sub_fasta(fasta, aln_path, rt_reg, "reg_RT")
    
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
                        uncover = fasta.count('-') + fasta.count('N')
                        line = pd.DataFrame(data = {"SAMPLE_No_NGS" : str(sample),
                                                    "Region_Protein" : gene,
                                                    'fasta': fasta,
                                                    "%coverage": round(100*((len(fasta)-uncover)/len(fasta)),2) if not uncover == len(fasta) else 0}
                                            , index=[0])
                        df = pd.concat([df,line], axis = 0)
        df.to_csv("QC/fasta.csv", index=False)
        format_df = pd.read_excel(forma,converters={'SAMPLE_No_NGS':str})
        mergi = pd.merge(format_df, df, on=["SAMPLE_No_NGS","Region_Protein"], how = "left")
        mergi.to_excel("QC/final_report.xlsx", index=False)
        

    #@override
    def cns(self, bam_path, cns_path, cns_x_path, min_depth_call, min_freq_thresh):
        
        vcf_path = bam_path.replace("BAM","VCF")
        for file in os.listdir(vcf_path):
            if file.endswith("csv"):
                vcf = pd.read_csv(vcf_path + file)
                sample = file.split(".csv")[0]
                vcf["base"] = vcf[["%A","%T","%C","%G"]].apply(lambda row: row[row > 20].nlargest(2).index.values, axis=1)
                #clean
                vcf["base"] = vcf["base"].astype(str).apply(lambda row: row.replace("%", "").replace("'","").replace("[","").replace("]",""))
                vcf.loc[vcf["base"] == "", "base"] = 'A C T G'
                deg_nuc = pd.read_csv(SCRIPT_PATH + "/refs/degenerate_nuc.csv", sep='\t')
                
                vcf = pd.merge(vcf, deg_nuc, on=["base"], how = "left")
                vcf.to_csv(vcf_path + sample + ".csv", index = False)
                
                cns = "".join(vcf["cns"].to_list())
               
                with open(cns_path +sample + ".fasta", 'w') as f:
                    f.write(">" + sample + '\n' + cns + '\n')
            
    
    def mafft(self, not_aligned, aligned):
        '''
        multi-fasta align.
        cat all consensus fasta sequences and run MAFFT. the implementation of MAFFT is in utils.

        '''
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS/*"), shell=True)
        mafft(self.reference, not_aligned, aligned)
    
    def qc_report(self, bam_path, depth_path, output_report):
        self.cut_genes("alignment/")
        self.excel_fasta("alignment/",self.metadata)
        fasta = pd.read_csv(output_report.replace("QC_report", "fasta.csv"), dtype={'SAMPLE_No_NGS': 'string'})
        
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
        writer.writerow(['sample', '%mapped','mapped_reads','total_reads','cov_bases','%coverage_RT','%coverage_PR',"%coverage_INT", 'mean_depth','max_depth','min_depth'])
        
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    sample = bam_file.split(".mapped")[0]
                    total_reads = pysam.AlignmentFile(bam_path+bam_file.split(".mapped")[0]+".bam").count(until_eof=True) #need to fix 
                    cover_rt = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "RT"), "%coverage"].reset_index(drop=True)
                    cover_rt = 0 if len(cover_rt) == 0 else cover_rt[0]
                    cover_pr = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "PR"), "%coverage"].reset_index(drop=True)
                    cover_pr = 0 if len(cover_pr) == 0 else cover_pr[0]
                    cover_int = fasta.loc[(fasta["SAMPLE_No_NGS"] == sample) & (fasta["Region_Protein"] == "Int"), "%coverage"].reset_index(drop=True)
                    cover_int = 0 if len(cover_int) == 0 else cover_int[0]
                    
                    
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