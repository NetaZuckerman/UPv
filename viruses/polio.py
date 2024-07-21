#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:16:08 2022

@author: hagar

polio pipeline considering all 3 lineages of polio, aiming to recognize the lineages exist in each sample. 
the assumtion here is that we get a mixture of polio viruses from several lineage in each sample. 
if you know the specific lineage in your sample, use the general pipeline.
"""


FILTER = "bwa mem -v1 -t %(threads)s %(reference)s %(r1)s %(r2)s | samtools view -@ 16 -b -f 2 > %(output_path)s%(sample)s.bam"
BAM2FQ = "samtools bam2fq -0n %(file)s.bam > %(file)s.fastq"

from utils import utils
from pipelines.generalPipeline import  general_pipe
import subprocess
import os

class polio(general_pipe):
    def __init__(self, reference, fastq, minion, threads):
        '''

        Parameters
        ----------
        reference : str
            path to the reference fasta file.
        fastq : str
            path to fastq folder
        minion : BOOL
            boolean to indicate if the the reads are minion based. defualt is illumina
        threads : int
            max number of threads for parts threads are available in this pipeline.

        Raises
        ------
        ValueError
            polio pipeline only works with illumina reads. for minion use PoP.

        Returns
        -------
        None.

        '''
        super().__init__(reference, fastq, minion, threads) 
        if self.minion:
            raise ValueError("polio pipeline only works with illumina reads. for minion use PoP.\nhttps://github.com/NetaZuckerman/PoP")
        utils.create_dirs([self.fastq+"polio_reads"])


    def filter_non_polio(self):
        '''
        filter outs reads that are not mapped to polio.
        
        Parameters
        ----------

        Returns
        -------
        None.

        '''
        for sample, r1 in self.sample_fq_dict.items():
            r2 = r1.replace("R1","R2")
            subprocess.call(FILTER % dict(threads = self.threads, reference=self.reference, r1=self.fastq + r1, r2=self.fastq + r2, output_path=self.fastq+"polio_reads/", sample=sample), shell=True)
            subprocess.call(BAM2FQ % dict(file=self.fastq+"polio_reads/" + sample), shell=True)
            os.remove(self.fastq+"polio_reads/" + sample + ".bam")
    
    def mapping(self):
        '''
        use general pipeline's mapping and split the result 
        so each sample will have a bam for each segment.
    
        Returns
        -------
        None.
    
        '''
        self.filter_non_polio()
        super().mapping()
        utils.split_bam("BAM/")

        
    
    def qc_report(self, bam_path, depth_path, output_report):

        super().results_report(bam_path, depth_path, output_report)
    