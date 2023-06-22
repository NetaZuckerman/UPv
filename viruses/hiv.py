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
from utils.utils import get_sequences
import os
SCRIPT_PATH = os.path.dirname(__file__)

class hiv(general_pipe):
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq)    
        
    def excel_fasta(self, cns_path):
        df = pd.DataFrame(columns=["sample", "gene", "fasta"])
        for file in listdir(cns_path):
            sample = file.split(".REF")[0]
            gene = file.split(".fa")[0].split("_")[-1]
            fasta = get_sequences(cns_path + file).popitem()[1]
            
            line = pd.DataFrame(data = {"sample" : sample,
                                 "gene" : gene,
                                 'fasta': fasta}, index=[0])
            df = pd.concat([df,line], axis = 0)
            df.to_csv("QC/fasta.csv", index=False)
                
        
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
            
            