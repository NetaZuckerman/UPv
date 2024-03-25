#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 11:52:43 2022

@author: hagar
"""

from pipelines.generalPipeline import general_pipe
from viruses.flu import flu
from viruses.flu_de_novo import flu_de_novo
from pipelines.deNovo import de_novo
from pipelines.multiRef import multi_ref
from viruses.polio import polio
from viruses.hiv import hiv
from viruses.hsv import hsv
from viruses.cmv import cmv

from utils import utils, parse_gb_file, parse_input

from mutations import signatures
import os
import sys


SCRIPT_PATH = os.path.dirname(__file__) +'/'

if __name__ == "__main__":
    
    DEBUG = False
    
    reference, fastq, threads, flu_flag, de_novo_flag, flu_de_novo_flag, polio_flag, cmv_flag, hiv_flag, \
        hsv_flag, minion_flag, gb_file, regions_file, mutations_flag, mini, \
            vcf, cns_min_freq_thresh, cns_min_depth_call, multi_ref_flag, drop_joint_reads = parse_input.parser()
    
    with open("log.txt",'w') as log:
        log.write("Working Directory: \n" + os.getcwd() + '\n\n')
        log.write("Command: \n" + ' '.join(sys.argv))
    
    dirs=['BAM','QC','CNS','CNS_' + cns_min_depth_call ,'QC/depth', 'alignment'] if not hiv_flag else ['BAM','QC','CNS','QC/depth', 'alignment']
    
    if hiv_flag:
        reference = SCRIPT_PATH + "viruses/refs/K03455.1_HIV.fasta"
        vcf = True
    elif not reference and not (hsv_flag or cmv_flag):
        raise ValueError("reference sequence is required.")
        
    if not mini:      
        utils.create_dirs(dirs) 
        
        if fastq and not mini:
            if not fastq.endswith("/"):
                fastq = fastq+"/"
        
            
            if flu_flag:
                pipe = flu(reference,fastq, minion_flag, threads)
            
            elif de_novo_flag:
                pipe = de_novo(reference,fastq, minion_flag, threads, de_novo_flag) 
            
            elif flu_de_novo_flag:
                pipe = flu_de_novo(reference, fastq, minion_flag, threads)
            
            elif polio_flag:
                pipe = polio(reference,fastq, minion_flag, threads)
            
            elif hiv_flag:
                pipe = hiv(reference,fastq, minion_flag, threads, hiv_flag)
                
            elif hsv_flag:
                pipe = hsv(reference,fastq, minion_flag, threads)
            
            elif cmv_flag:
                pipe = cmv(reference,fastq, minion_flag, threads)

            elif multi_ref_flag:
                pipe = multi_ref(reference,fastq, minion_flag, threads, drop_joint_reads)
    
            else:
                pipe = general_pipe(reference,fastq, minion_flag, threads)
            
           
            # mapping 
            pipe.mapping()
            
        
            if vcf:
                utils.create_dirs(["VCF"])
                pipe.variant_calling('BAM/', 'VCF/')
                
            pipe.cns("BAM/","CNS/","CNS_" + cns_min_depth_call + '/', cns_min_depth_call, cns_min_freq_thresh) #temp comment
            
            pipe.depth("BAM/","QC/depth/")
            
            pipe.mafft("alignment/all_not_aligned.fasta", "alignment/all_aligned.fasta")
            
            pipe.qc_report("BAM/", "QC/depth/", 'QC/QC_report') #temp comment
            
    else: #yes mini
        if not os.path.exists("alignment"):    
            utils.create_dirs(["alignment"])
        utils.mafft(reference, fastq, "alignment/all_aligned.fasta")
        
    if gb_file:
        parse_gb_file.parse(gb_file)
        regions_file = gb_file.replace(".gb", "_regions.csv")

    if mutations_flag or mini and not cmv_flag:
        utils.create_dirs(["reports"])
        signatures.run("alignment/all_aligned.fasta", regions_file, "reports/mutations.xlsx")
    
    