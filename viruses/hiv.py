#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 06:40:21 2023

@author: hagar

This script is designed specifically for our HIV department (Orna).
The illumina reads are mapped to the whole genome reference.
The genes: RT(reverese transcriptasem), PR(protease), int(integrase) are cut.
consensus sequence of each sample is determined by Orna's requermients and explained in cns().
The final report is a metadata file that HIV department suppose to provide us, merged with the consensus sequences. 
"""

import pandas as pd
from pipelines.generalPipeline import general_pipe, ALL_NOT_ALIGNED
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
BOWTIE_INDEX = "bowtie2-build %(reference)s.fasta %(reference)s"
BOWTIE2 = "bowtie2 -x %(reference)s -1 %(r1)s -2 %(r2)s --threads %(threads)s  --very-sensitive-local | samtools view -@ %(threads)s -b - > %(output_path)s%(sample)s.bam"
#*gene_regions*#
pr_reg = (2252,2550)
rt_reg = (2661,3294)
int_reg = (4230, 5094)

class hiv(general_pipe):
    def __init__(self, reference, fastq,minion, threads, metadata, sensitive):
        '''
        

        Parameters
        ----------
        reference : str
            path to the reference fasta file.
            no need to provide in HIV runs.
            using K03455.1 as reference. it is defined in the main script.
        fastq : str
            path to fastq folder
        minion : BOOL
            boolean to indicate if the the reads are minion based. defualt is illumina
        threads : int
            max number of threads for parts threads are available in this pipeline.
        metadata : str
            path to xlsx file HIV department should send as.

        Returns
        -------
        None.

        '''
        super().__init__(reference, fastq, minion, threads)    
        self.metadata = metadata
        self.sensitive = sensitive
 
    
    def cut_genes(self, aln_path):
        '''
        cut the genes from the all_aligned file.

        Parameters
        ----------
        aln_path : str
            path to alignment fasta file.

        Returns
        -------
        None.

        '''
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
        
    def mapping(self):
        
        if self.sensitive:
            subprocess.call(BOWTIE_INDEX % dict(reference=self.reference.split(".fa")[0]), shell=True)
            for sample, fq in self.sample_fq_dict.items():
                r1 = fq
                r2 = r1.replace("R1","R2")
                subprocess.call(BOWTIE2 % dict(threads = self.threads, reference=self.reference.split(".fa")[0],r1=self.fastq + r1, r2=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
                super().process_bam(sample)
        else:
            super().mapping()
            
    # def cns(self, bam_path, cns_path, cns_x_path, min_depth_call, min_freq_thresh):
    #     '''
    #     @override        
        
    #     generate consesus sequence from VCF files.
    #     rules:
            
        
    #     Parameters
    #     ----------

    #     Returns
    #     -------
    #     None.

    #     '''        
    #     #override parameters 
    #     min_freq_thresh = 10 #the minimal %depth of a base to be considered
    #     min_depth_call = 0 # the minimal read depth in each position
        
    #     vcf_path = bam_path.replace("BAM","VCF")
    #     for file in os.listdir(vcf_path):
    #         if file.endswith("csv"):
    #             vcf = pd.read_csv(vcf_path + file)
    #             sample = file.split(".csv")[0]
                
    #             #to determine consensus choose the 2 bases with depth > 5% 
    #             #nlargest(2): the 2 largest depths 
    #             #save the chosen bases in column "base"
    #             vcf["base"] = vcf[["%A","%T","%C","%G"]].apply(lambda row: row[row > int(min_freq_thresh)].nlargest(2).index.values, axis=1)
                
    #             #clean unwanted characters
    #             vcf["base"] = vcf["base"].astype(str).apply(lambda row: row.replace("%", "").replace("'","").replace("[","").replace("]",""))
                
    #             #if "base" is empty it means that no base had depth > 5%. set these cases to 'A C T G'
    #             vcf.loc[vcf["base"] == "", "base"] = 'A C T G'
                
    #             #if depth <= min_depth_call set it as 'A C T G'
    #             vcf.loc[vcf["depth"] <= int(min_depth_call), "base"] = 'A C T G'
                
    #             #this table contains the translation of "base" column to consensus nucleotide.
    #             deg_nuc = pd.read_csv(SCRIPT_PATH + "/refs/degenerate_nuc.csv", sep='\t')
                
    #             #add the translation to consensus nucleotide
    #             vcf = pd.merge(vcf, deg_nuc, on=["base"], how = "left")
    #             vcf.to_csv(vcf_path + sample + ".csv", index = False)
                
    #             #convert consensus to string
    #             cns = "".join(vcf["cns"].to_list())
               
    #             #save consensus
    #             with open(cns_path +sample + ".fasta", 'w') as f:
    #                 f.write(">" + sample + '\n' + cns + '\n')
            
    
    def mafft(self, not_aligned, aligned):
        '''
        multi-fasta align.
        cat all consensus fasta sequences and run MAFFT. the implementation of MAFFT is in utils.

        '''
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS/*"), shell=True)
        mafft(self.reference, not_aligned, aligned)
    
    def qc_report(self, bam_path, depth_path, output_report):
        '''
        generate qc report considering each gene information.
        generate "final_report" using excel_fasta()

        Parameters
        ----------
        bam_path : str
            path to bam folder.
        depth_path : TYPE
           path to depth folder
        output_report : TYPE
            path to qc report file

        Returns
        -------
        None.

        '''
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