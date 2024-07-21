#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 10:15:00 2022

@author: hagar

analyse known subtype of influenza. 
the provided reference should be a multi-fasta of the 8 segments of flu.


"""
from pipelines.generalPipeline import general_pipe, ALL_NOT_ALIGNED
import os
import subprocess
from utils import utils

#bash commands
SPLIT = "bamtools split -in %(bam)s -reference"
CAT_SAMPLE = "cat CNS_5/%(sample)s.REF* > CNS_5/per_sample/%(sample)s.fasta"
CAT_NO_HEADER = "awk 'FNR>1' %(files)s* >> %(bigfile)s"

class flu (general_pipe):

    def __init__(self, reference, fastq, minion, threads):
        super().__init__(reference, fastq, minion, threads)    
   
    def mapping(self):
        '''
        use general pipeline's mapping and split the result 
        so each sample will have a bam for each segment.

        Returns
        -------
        None.

        '''
        super().mapping()
        utils.split_bam("BAM/")
        
    def concat_samples(self):
        '''
        cns() generated one seperated fasta files for each sample segment
        generate a multi-fasta file for each sample with all segments.
        save in cns_path/per_sample/ folder.
        save cns_path/per_sample/cat/ concatenated version (one header) for later mafft.

        Returns
        -------
        None.

        '''
        os.makedirs("CNS_5/per_sample")
        os.makedirs("CNS_5/per_sample/cat/")
        skip_files=["R2", "Undetermined", "unpaired", "singletons"]
        for sample in os.listdir(self.fastq):
            if any(skip_file in sample for skip_file in skip_files):
                continue
            sample = sample.split("_")[0] #sample short name
            subprocess.call(CAT_SAMPLE % dict(sample=sample), shell=True)
            #create file with no headers (for mafft)
            cat_file_name = "CNS_5/per_sample/cat/" + sample + ".fa"
            file = open(cat_file_name, 'w')
            file.write(">"+sample+'\n')
            file.close()
            subprocess.call(CAT_NO_HEADER % dict(files="CNS_5/"+sample,bigfile=cat_file_name), shell=True)
        
    def concat_segments(self):
        '''
        generate a multi-fasta file for each segment with all samples.

        Returns
        -------
        None.

        '''
        os.makedirs("CNS_5/per_gene") 
        #get segments list 
        ref= open(self.reference,"r")
        segments = []
        for line in ref:
            if line.startswith(">"):
                segments.append(line[1:-1])
        for segment in segments:
            segment_file = open("CNS_5/per_gene/" + segment + ".fa","a+")
            for file_name in os.listdir("CNS_5/"):        
                 if segment in file_name:
                     file = open("CNS_5/" + file_name, 'r')
                     segment_file.write(file.read())
                     file.close()
        ref.close()
        segment_file.close()
        
            
    def mafft(self, not_aligned, aligned):
        '''
        align all concatnated samples to concatenated reference.
        
        Parameters
        ----------
        not_aligned : str
            path to not aligned fasta file.
        aligned : str
            path to aligned fasta file.
            
        Returns
        -------
        None.

        '''
        self.concat_samples()
        self.concat_segments()
        
        #all samples without segment and not aligned
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS_5/per_sample/cat/*"), shell=True)  
        #remove segments headers from reference and save it
        if not os.path.isfile("reference_for_mafft.fa"):
            reference_for_mafft = open("reference_for_mafft.fa", 'a+')
            reference = open(self.reference,'r')
            reference_for_mafft.write(">reference \n")
            lines=reference.readlines()
            lines = [s.rstrip('\n') for s in lines]
            for line in lines:
                if not line.startswith(">"):
                    reference_for_mafft.write(line)
            reference_for_mafft.close()

        #run mafft
        utils.mafft("reference_for_mafft.fa", "alignment/all_not_aligned.fasta", "alignment/all_aligned.fasta")
        os.remove("reference_for_mafft.fa")