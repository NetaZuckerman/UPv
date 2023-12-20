#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 06:43:36 2023

@author: hagar
"""

from pipelines.generalPipeline import general_pipe
import os
import pandas as pd
import subprocess

UNIQ_READS = "samtools view -bq 1 %(bam_file)s > %(ouput)s"
SPLIT = "bamtools split -in %(bam)s -reference"


class multi_ref (general_pipe):

    def __init__(self, reference, fastq, threads, drop_joint_reads= False):
        super().__init__(reference, fastq, threads)    
        self.drop_joint_reads = drop_joint_reads
   
    def mapping(self):
        super().mapping()
        
        bam_path = "BAM/"
        #remove reads mapped to multiple reference 
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file and "bai" not in bam_file:
                if self.drop_joint_reads:
                    uniq_file = bam_file.replace("bam", "uniq.bam")
                    subprocess.call(UNIQ_READS % dict(bam_file = bam_path + bam_file, ouput = bam_path + uniq_file), shell=True)
                    subprocess.call(SPLIT % dict(bam=bam_path + uniq_file), shell=True)#split bams by reference
                    os.remove(bam_path + uniq_file)
                else:
                    subprocess.call(SPLIT % dict(bam=bam_path + bam_file), shell=True)#split bams by reference
                os.remove(bam_path + bam_file)
                
        
        
    #skip mafft in multi reference pipeline
    def mafft(self, not_aligned, aligned):
        return
    
    
    def qc_report(self, bam_path, depth_path, output_report):
        super().qc_report("BAM/", "QC/depth/", 'QC/QC_report')
        qc = pd.read_csv("QC/QC_report.csv")
        qc[['sample', 'reference']] = qc['sample'].str.split('.REF_', expand=True)
        #reorder columns
        qc = qc[['sample', 'reference','mapped%','mapped_reads','total_reads','cov_bases','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth', 'chimeric_read_count']]
        qc.to_csv("QC/QC_report.csv", index=False)