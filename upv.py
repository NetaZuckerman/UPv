#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:52:43 2022

@author: hagar
"""

from pipelines.generalPipeline import general_pipe
from viruses.flu import flu
from pipelines.deNovo import de_novo
from pipelines.multiRef import multi_ref
from viruses.polio import polio
from viruses.hiv import hiv
from viruses.hsv import hsv
from utils import utils, parse_gb_file, parse_input
from mutations import signatures
import os


SCRIPT_PATH = os.path.dirname(__file__) +'/'

if __name__ == "__main__":
    
    DEBUG = False
    
    reference, fastq, threads, flu_flag, de_novo_flag, polio_flag, cmv_flag, hiv_flag, \
        hsv_flag, gb_file, regions_file, mutations_flag, mini, \
            vcf, cnsThresh, multi_ref_flag, drop_joint_reads = parse_input.parser()
    
    dirs=['BAM','QC','CNS','CNS_5','QC/depth', 'alignment'] if not hiv_flag else ['BAM','QC','CNS','QC/depth', 'alignment']
    
    if hiv_flag:
        reference = SCRIPT_PATH + "viruses/refs/K03455.1_HIV.fasta"
        vcf = True
    elif not reference:
        raise ValueError("reference sequence is required.")
        
    if not mini:      
        utils.create_dirs(dirs) 
        
        if fastq and reference and not mini:
            if not fastq.endswith("/"):
                fastq = fastq+"/"
        
            
            if flu_flag:
                pipe = flu(reference,fastq, threads)
            
            elif de_novo_flag:
                pipe = de_novo(reference,fastq, threads, de_novo_flag) 

            elif polio_flag:
                pipe = polio(reference,fastq, threads)
            
            elif hiv_flag:
                pipe = hiv(reference,fastq, threads, hiv_flag)
                
            elif hsv_flag:
                pipe = hsv(reference,fastq, threads)
            
            elif multi_ref_flag:
                pipe = multi_ref(reference,fastq, threads, drop_joint_reads)
    
            else:
                pipe = general_pipe(reference,fastq, threads)
            
           
            #mapping 
            pipe.mapping()
            
        
            if vcf:
                utils.create_dirs(["VCF"])
                pipe.variant_calling()
                
            pipe.cns_depth("BAM/","QC/depth/","CNS/","CNS_5/", cnsThresh) #temp comment
            
            pipe.mafft("alignment/all_not_aligned.fasta", "alignment/all_aligned.fasta")
       
            
            pipe.qc_report("BAM/", "QC/depth/", 'QC/QC_report') #temp comment
            

    else: #yes mini
        utils.create_dirs(["alignment"])
        utils.mafft(reference, fastq, "alignment/all_aligned.fasta")
        
    if gb_file:
        parse_gb_file.parse(gb_file)
        regions_file = gb_file.replace(".gb", "_regions.csv")

    if mutations_flag or mini:
        if not gb_file and not regions_file:
            raise ValueError("gene bank file or regions file is required.")
        utils.create_dirs(["reports"])
        signatures.run("alignment/all_aligned.fasta", regions_file, "reports/mutations.xlsx")
    
    