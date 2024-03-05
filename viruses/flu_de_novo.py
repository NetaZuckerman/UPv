#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 10:21:47 2024

@author: hagar
"""

from pipelines.generalPipeline import general_pipe, INDEX, BWM_MEM, FILTER_BAM, SORT,CHIMER
from pipelines.deNovo import de_novo
from utils.utils import create_dirs, split_bam
import pandas as pd
from Bio import SeqIO
import os
import subprocess

class flu_de_novo(general_pipe):
      
    def __init__(self, reference, fastq, minion, threads):
        super().__init__(reference, fastq, minion, threads) 
        if self.minion:
            raise ValueError("de-novo pipeline only works with illumina reads.")
        
    
    def choose_references(self):
        create_dirs(["fasta/","fasta/selected_contigs","fasta/all_contigs","fasta/selected_references"])

        for file in(os.listdir("blast/")):
            if file.endswith("csv") and "filter" not in file:
                sample = file.split(".csv")[0]
                blast = pd.read_csv("blast/" + file)
                blast = blast[blast['subject_title'].str.contains('Influenza A')]
                  
                pb2 =  blast[blast['subject_title'].str.contains('PB2|segment 1')]
                pb1 =  blast[blast['subject_title'].str.contains('PB1|segment 2')]
                pa =  blast[blast['subject_title'].str.contains('PA|segment 3')]
                ha =  blast[blast['subject_title'].str.contains('HA|segment 4')]
                np =  blast[blast['subject_title'].str.contains('NP|segment 5')]
                na =  blast[blast['subject_title'].str.contains('NA|segment 6')]
                m =  blast[blast['subject_title'].str.contains('M1|M2|segment 7')]
                ns =  blast[blast['subject_title'].str.contains('NS1|NS2|segment 8')]
                
                segments = {"pb2": pb2, "pb1": pb1, "pa": pa, "ha": ha, "np": np, "na": na, "m": m, "ns": ns}
                
                for segment, df in segments.items():
                    reference = df[df['raw_score'] == df['raw_score'].max()]
                    if len(reference)>0:
                        reference_id = reference["subject_seq_id"].iloc[0]
                        reference_name = reference["subject_title"].iloc[0]
                        if '|' in reference_id:
                            reference_id = reference_id.split('|')[1]
                        for record in SeqIO.parse(self.reference, "fasta"):
                            # Check if the current record's header matches the desired header
                            if reference_id in record.description:
                                with open("fasta/selected_references/" + sample + ".fasta", 'a') as f:
                                    f.write(">" + segment + "_" + reference_id + '_' + reference_name + '\n' + str(record.seq) + '\n')

        
    def mapping(self):
        # de_novo.run_spades(self)
        # de_novo.run_blast(self)
        self.choose_references()
        for sample, r1 in self.sample_fq_dict:
            r2 = r1.replace("R1","R2")
            fasta = sample + ".fasta"
            if os.path.exists("fasta/selected_references/" + fasta):
                subprocess.call(INDEX % dict(reference="fasta/selected_references/" + fasta), shell=True)
                subprocess.call(BWM_MEM % dict(reference="fasta/selected_references/" + fasta ,r1=self.fastq+r1, r2=self.fastq+r2, sample=sample, output_path="BAM/",threads = self.threads), shell=True) #map to reference
                subprocess.call(FILTER_BAM % dict(sample=sample,filter_out_code = 4, output_path="BAM/",threads = self.threads), shell=True) #keep mapped reads
                subprocess.call(SORT % dict(sample=sample, output_path="BAM/",threads = self.threads), shell=True)         
                subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/"), shell=True)
                split_bam("BAM/")
        
    #override
    def mafft(self, not_aligned, aligned):
        return
        