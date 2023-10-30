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
BWM_MEM = "bwa mem -v1 -t %(threads)s %(reference)s %(r1)s %(r2)s | samtools view -@ %(threads)s -b - > %(output_path)s%(sample)s.bam"
FILTER_BAM = "samtools view -@ %(threads)s -b -F %(filter_out_code)s %(output_path)s%(sample)s.bam > %(output_path)s%(sample)s.mapped.bam"
SORT = "samtools sort -@ %(threads)s %(output_path)s%(sample)s.mapped.bam -o %(output_path)s%(sample)s.mapped.sorted.bam"
CHIMER = "samtools view %(output_path)s%(sample)s.bam |  grep 'SA:' > %(output_path)s%(sample)s.chimeric_reads.txt"
DEPTH = "samtools depth -a %(bam_path)s%(bam_file)s > %(depth_path)s%(sample)s.txt"
CNS = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t %(cnsThresh)s -m 1 -p %(cns_path)s%(sample)s.fa"
CNS5 = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t %(cnsThresh)s -m 5 -p %(cns5_path)s%(sample)s.fa"
#variants
VCF = "bcftools mpileup -d 8000 --skip-indels -f %(ref)s %(bam)s -a AD,DP | sed '/^#/d' > %(vcf)s"


#MAFFT commands
ALL_NOT_ALIGNED =   "cat %(dir)s > alignment/all_not_aligned.fasta"

#report commands
BREADTH_CNS5 = "$(cut -f3 QC/depth/%(sample)s.txt | awk '$1>5{c++} END{print c+0}')"
#directories for output

            
class general_pipe():
    
    
    def __init__(self, reference, fastq, threads):
        self.reference = reference
        self.fastq = fastq
        self.r1r2_list = utils.get_r1r2_list(self.fastq)
        self.threads = threads
        subprocess.call(INDEX % dict(reference=self.reference), shell=True)


    def bam(self):
        '''
        generate bam file from paired-end fastq (R1,R2).
        generate filtered file with mapped reads.
        generate sorted file
        
        Parameters
        ----------
        sample_r1_r2 : list of lists
            [sample, r1 fastq path, r2 fastq path].

        '''
        filter_out_code = 2052
         
        for sample, r1, r2 in self.r1r2_list:
            subprocess.call(BWM_MEM % dict(threads = self.threads, reference=self.reference,r1=self.fastq + r1, r2=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
            subprocess.call(FILTER_BAM % dict(threads = self.threads, sample=sample, filter_out_code = filter_out_code, output_path="BAM/"), shell=True) #keep reads according to the filter code: https://broadinstitute.github.io/picard/explain-flags.html
            subprocess.call(SORT % dict(threads = self.threads, sample=sample, output_path="BAM/"), shell=True)         
            #chimeric reads: find reads that where mapped to more then one region in the genome
            subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/"), shell=True)
        
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
            if "sorted" in bam:
                sample = bam.split(".bam")[0].replace(".mapped.sorted", "")
                vcf_file = vcf_path + sample + ".vcf"
                subprocess.call(VCF % dict(bam=bam_path + bam, ref=self.reference, vcf= vcf_file ), shell=True) 
                vcf = pd.read_csv(vcf_file, sep='\t', names=["ref_id",	"pos",	"id",	"ref",	"alt",	"qual"	,"filter",	"info"	,"format", "details"])
                #process alt nucleotides accurance
                alts = vcf["alt"].str.split(",",expand=True)
                # add column if not exist
                if 2 not in alts.columns:
                    alts[2] = None
                cols = ["alt" + str(x) for x in alts.columns]
                alts.set_axis(cols, axis=1,inplace=True)
                vcf = pd.concat([vcf, alts], axis=1)
                
                #process each nucleotide depth
                depths = vcf["details"].str.split(":",expand=True)[2].str.split(",",expand=True)
                # add column if not exist
                if 3 not in depths.columns:
                    depths[3] = None
                
                
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
                utils.fix_cns_header("CNS/")
                utils.fix_cns_header("CNS_5/")
    def mafft(self, not_aligned, aligned):
        '''
        multi-fasta align.
        cat all consensus fasta sequences and run MAFFT. the implementation of MAFFT is in utils.

        '''
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS_5/*"), shell=True)
        utils.mafft(self.reference, not_aligned, aligned)
    
    #write report.csv - mapping analysis
    #de novo flag was created to generate different output for de novo analysis, used in denovo class and polio class.

    def depth(self, file, start=0, end=-1):
        depths = [int(x.split('\t')[2]) for x in open(file).readlines()]
        
        #set
        if end == -1:
            genome_size = len(depths)   
        else:
            depths = depths[start:end]
            genome_size = end-start
        #calculate depths
        depths = [i for i in depths if i != 0]
        if not depths:
            return '','','','',''
        max_depth = max(depths)
        min_depth = min(depths)
        mean_depth = str(round(mean(depths),3))
        coverage = round(len(depths) / genome_size * 100,3)  if genome_size else ''
        
        #cns 5: position with depth > 5
        breadth_cns5 = len([i for i in depths if i > 5])
        cns5_cover = round(breadth_cns5 / genome_size * 100,3)  if genome_size else ''
        
        return max_depth, min_depth,mean_depth, coverage, cns5_cover
    
    def chimer_count(self,file_name):
        #count chimeric reads
        with open(file_name) as f:
            chimer_count = len(f.readlines())
        return chimer_count
    
    def general_qc(self, bam_file):
        total_reads = pysam.AlignmentFile(bam_file.split(".mapped")[0]+".bam").count(until_eof=True) 
        coverage_stats = pysam.coverage(bam_file).split("\t")
        mapped_reads = int(coverage_stats[11])
        mapped_percentage = round(mapped_reads/total_reads*100,4) if total_reads else ''
        cov_bases =  int(coverage_stats[12])
        chimer_count = self.chimer_count(bam_file.split(".mapped")[0].split(".REF")[0] + ".chimeric_reads.txt")
        
        return total_reads, mapped_reads, mapped_percentage, cov_bases, chimer_count
        
    def qc_report(self, bam_path, depth_path, output_report):
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
        writer.writerow(['sample', 'mapped%','mapped_reads','total_reads','cov_bases','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth', 'chimeric_read_count'])
        
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    
                    #general qc
                    total_reads, mapped_reads, mapped_percentage, cov_bases, chimer_count = self.general_qc(bam_path+bam_file)
                    
                    #depth 
                    max_depth, min_depth, mean_depth, cover, cns5_cover = self.depth(depth_path+sample+".txt")
                                        
                    writer.writerow([sample, mapped_percentage, mapped_reads, total_reads, cov_bases, cover, cns5_cover, mean_depth, max_depth, min_depth, chimer_count])
                
        f.close()
    def mafft(self, not_aligned, aligned):
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS_5/*"), shell=True)  
        utils.mafft(self.reference, not_aligned, aligned)