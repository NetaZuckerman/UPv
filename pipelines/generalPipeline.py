#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:37:02 2022

@author: hagar

General pipeline is the base class of the upv pipeline.
"""

import subprocess
import os
import pysam
from statistics import mean
import csv
from utils import utils
import pandas as pd


#alignment conmmands
INDEX = "bwa index %(reference)s"
#BWM_MEM = "bwa mem -v1 -t 16 %(reference)s %(r1)s | samtools view -@ 16 -b - > %(output_path)s%(sample)s.bam"
BWM_MEM = "bwa mem -v1 -t 16 %(reference)s %(r1)s %(r2)s | samtools view -@ 16 -b - > %(output_path)s%(sample)s.bam"
MAPPED_BAM = "samtools view -@ 16 -b -F 4 %(output_path)s%(sample)s.bam > %(output_path)s%(sample)s.mapped.bam"
SORT = "samtools sort -@ 14 %(output_path)s%(sample)s.mapped.bam -o %(output_path)s%(sample)s.mapped.sorted.bam"
DEPTH = "samtools depth -a %(bam_path)s%(bam_file)s > %(depth_path)s%(sample)s.txt"
CNS = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t %(cnsThresh)s -m 1 -p %(cns_path)s%(sample)s.fa"
CNS5 = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t %(cnsThresh)s -m 5 -p %(cns5_path)s%(sample)s.fa"
SAMTOOLS_INDEX = "samtools index -@ 16 %(bam_path)s%(bam_file)s"
#variants
VCF = "bcftools mpileup -d 8000 --skip-indels -f %(ref)s %(bam)s -a AD,DP | sed '/^#/d' > %(vcf)s"


#MAFFT commands
ALL_NOT_ALIGNED =   "cat %(dir)s > alignment/all_not_aligned.fasta"

#report commands
BREADTH_CNS5 = "$(cut -f3 QC/depth/%(sample)s.txt | awk '$1>5{c++} END{print c+0}')"
#directories for output

            
class general_pipe():
    
    
    def __init__(self, reference, fastq):
        self.reference = reference
        self.fastq = fastq
        self.r1r2_list = utils.get_r1r2_list(self.fastq)
        subprocess.call(INDEX % dict(reference=self.reference), shell=True)


    def bam(self, sample_r1_r2):
        '''
        generate bam file from paired-end fastq (R1,R2).
        generate filtered file with mapped reads.
        generate sorted file
        
        Parameters
        ----------
        sample_r1_r2 : list of lists
            [sample, r1 fastq path, r2 fastq path].

        '''
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2]            
        #subprocess.call(BWM_MEM % dict(reference=self.reference,r1=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
        subprocess.call(BWM_MEM % dict(reference=self.reference,r1=self.fastq + r1, r2=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
        subprocess.call(MAPPED_BAM % dict(sample=sample, output_path="BAM/"), shell=True) #keep mapped reads
        subprocess.call(SORT % dict(sample=sample, output_path="BAM/"), shell=True)         
    
    def variant_calling(self, bam_path = "BAM/", vcf_path= "VCF/"):
        '''

        Parameters
        ----------
        bam_path : str, optional
            The default is "BAM/".
        vcf_path : str, optional
            The default is "VCF/".

        '''
        
        
        for bam in os.listdir(bam_path):
            if "REF" in bam:
                sample = bam.split(".bam")[0].replace(".mapped.sorted", "")
                vcf_file = vcf_path + sample + ".vcf"
                subprocess.call(VCF % dict(bam=bam_path + bam, ref=self.reference, vcf= vcf_file ), shell=True) 
                vcf = pd.read_csv(vcf_file, sep='\t', names=["ref_id",	"pos",	"id",	"ref",	"alt",	"qual"	,"filter",	"info"	,"format", "details"])
                #process alt nucleotides accurance
                alts = vcf["alt"].str.split(",",expand=True)
                cols = ["alt" + str(x) for x in alts.columns]
                alts.set_axis(cols, axis=1,inplace=True)
                vcf = pd.concat([vcf, alts], axis=1)
                
                #process each nucleotide depth
                depths = vcf["details"].str.split(":",expand=True)[2].str.split(",",expand=True)
                cols = ["depth" + str(x) for x in depths.columns]
                depths.set_axis(cols, axis=1,inplace=True)
                vcf = pd.concat([vcf, depths], axis=1)

                #get A G C T counts
                final_df = pd.DataFrame(columns=["position", "ref", "A",  "G", "C", "T"])
                for i, row in vcf.iterrows():
                    counts = {"position":row["pos"],
                     "ref": row["ref"],
                        row["ref"]:row["depth0"],
                     row["alt0"]: row["depth1"],
                     row["alt1"]: row["depth2"],
                     row["alt2"]: row["depth3"]}
                    if "<*>" in counts.keys():
                        counts.pop("<*>")
                    if None in counts.keys():
                        counts.pop(None)
                    
                    counts = pd.DataFrame(data = counts, index=[0])
                    final_df = pd.concat([final_df,counts], axis = 0)
                    final_df = final_df.fillna(0)
                    
                #calc the frequencies
                final_df[["A","G","C","T"]] = final_df[["A","G","C","T"]].astype(int)
                final_df["depth"] = final_df["A"] + final_df["G"] +  final_df["C"] + final_df["T"]
                final_df["%A"] = round(final_df["A"]/final_df["depth"]*100)
                final_df["%G"] = round(final_df["G"]/final_df["depth"]*100)
                final_df["%C"] = round(final_df["C"]/final_df["depth"]*100)
                final_df["%T"] = round(final_df["T"]/final_df["depth"]*100)
                
                final_df.to_csv(vcf_path + sample + ".csv", index = False)
                
    #find mapping depth and consensus sequence 
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path, cnsThresh):
        '''
        Generate consensus sequences and calculate aligning depths from a bam file.
        '''
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file and "bai" not in bam_file:
                sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                subprocess.call(DEPTH % dict(bam_path=bam_path, bam_file=bam_file, depth_path=depth_path, sample=sample), shell=True) 
                #consensus
                subprocess.call(CNS % dict(bam_path=bam_path, bam_file=bam_file, cns_path=cns_path, sample=sample, cnsThresh=cnsThresh), shell=True) 
                subprocess.call(CNS5 % dict(bam_path=bam_path, bam_file=bam_file, cns5_path=cns5_path, sample=sample, cnsThresh=cnsThresh), shell=True)
                #remove qual files
                os.remove(cns_path+sample+".qual.txt")
                os.remove(cns5_path+sample+".qual.txt")

    def mafft(self):
        '''
        multi-fasta align.
        cat all consensus fasta sequences and run MAFFT. the implementation of MAFFT is in utils.

        '''
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS_5/*"), shell=True)
        utils.mafft(self.reference)
    
    #write report.csv - mapping analysis
    #de novo flag was created to generate different output for de novo analysis, used in denovo class and polio class.
    def results_report(self, bam_path, depth_path, output_report):
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
        writer.writerow(['sample', 'mapped%','mapped_reads','total_reads','cov_bases','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth'])
        
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    # subprocess.call(SAMTOOLS_INDEX % dict(bam_path=bam_path, bam_file=bam_file.split(".mapped")[0]+".bam"), shell=True) #to fix
                    total_reads = pysam.AlignmentFile(bam_path+bam_file.split(".mapped")[0]+".bam").count(until_eof=True) #need to fix 
                    coverage_stats = pysam.coverage(bam_path+bam_file).split("\t")
                    mapped_reads = int(coverage_stats[11])
                    mapped_percentage = round(mapped_reads/total_reads*100,4) if total_reads else ''
                    cov_bases =  int(coverage_stats[12])
                    coverage = float(coverage_stats[13])
            
                    #depth 
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    depths = [int(x.split('\t')[2]) for x in open(depth_path+sample+".txt").readlines()]
                    depths = [i for i in depths if i != 0]
                    mean_depth = str(round(mean(depths),3)) if depths else ''
                    min_depth = min(depths) if depths else ''
                    max_depth = max(depths) if depths else ''
                    
                    breadth_cns5 = len([i for i in depths if i > 5])
                    genome_size = sum(1 for line in open(depth_path+sample+".txt"))
                    coverage_cns5 = round(breadth_cns5 / genome_size * 100,3)  if genome_size else ''
                                        
                    writer.writerow([sample, mapped_percentage, mapped_reads, total_reads, cov_bases, coverage, coverage_cns5, mean_depth, max_depth, min_depth])
                
        f.close()
        

