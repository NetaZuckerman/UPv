#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:26:50 2022

@author: hagar
"""
import os
import subprocess
from utils import utils
import sys

#fastq_path = "/mnt/data/projects/Merav/VDPV3_2022/fastq/raw/"
#output_path = "/home/hagar/polio_tests/filter_s1_map_to_A20/filtered_reads/"
polio_refs = "/home/hagar/refs/polio/polio.fasta"
fastq_path = sys.argv[1]
output_path = sys.argv[2]

FILTER = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b -F 2 - | samtools fastq > %(output_path)s%(sample)s.fastq"
FILTER2 = "bwa mem  -v1 -t 32 %(reference)s %(r1)s %(r2)s| samtools view -@ 32 -b -f 4 -f 8 - | samtools bam2fq > %(output_path)s%(sample)s.fastq"
FILTER3 = "bwa mem -v1 -t 32 %(reference)s %(r1)s | samtools view -@ 32 -b -f 4 -  | samtools bam2fq > %(output_path)s%(sample)s.fastq"
RUN_SPADES = "spades --12 %(fastq)s -o %(output_path)s --rna"
def filter_with_singeltons():
    for fastq in os.listdir(fastq_path):
        if "R1" in fastq:
            r1 = fastq_path+fastq
            r2 = r1.replace("R1","R2")
            sample = fastq.split("_")[0]
            subprocess.call(FILTER % dict( reference=polio_refs, r1=r1, r2=r2, output_path=output_path, sample=sample), shell=True)
        
#run spades on paired end fastq files. the list must contain sample, r1 r2 file names 
def run_spades(fastq):
        sample = fastq.split(".fastq")[0].split("_")[0]
        utils.create_dirs(["spades/spades_results/"+sample])
        subprocess.call(RUN_SPADES % dict( fastq=fastq_path+fastq, output_path="/mnt/data/projects/Merav/girl2/no_polio_de_novo_analysis/spades/spades_results/"+ sample + "/"), shell=True)
         
def filter_unmapped_paired(fastq):
    sample = fastq.split(".fastq")[0].split("_")[0]
    r1 = fastq_path+fastq
    r2 = r1.replace("R1","R2")
    subprocess.call(FILTER2 % dict( reference=polio_refs, r1=r1, r2=r2, output_path=output_path, sample=sample), shell=True)
        
def filter_only_singletons():
    for fastq in os.listdir(fastq_path):
        if "R1" in fastq:
            r1 = fastq_path+fastq
            r2 = r1.replace("R1","R2")
            sample = fastq.split("_")[0]
            subprocess.call(FILTER3 % dict( reference=polio_refs, r1=r1, output_path=output_path, sample=sample), shell=True)

fastq_list=[]
for fastq in os.listdir(fastq_path):
    if "R1" in fastq:
        fastq_list.append(fastq)
#utils.run_mp(5,run_spades,fastq_list)
utils.run_mp(2,filter_unmapped_paired,fastq_list)
#filter_with_singeltons()
#filter_only_singletons()