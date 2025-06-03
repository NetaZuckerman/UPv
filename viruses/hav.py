#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 09:16:26 2025

@author: hagar

*no minion option*

"""

import os
from pipelines.generalPipeline import general_pipe
from pipelines.multiRef import SPLIT
import subprocess

SCRIPT_PATH = os.path.dirname(__file__) +'/'

class hav(general_pipe):
    def __init__(self, reference, fastq, minion, threads):
        '''

        Parameters
        ----------
        reference : str
            path to the reference fasta file. no need to provide reference.
        fastq : str
            path to fastq folder
        minion : BOOL
            boolean to indicate if the the reads are minion based. defualt is illumina
        threads : int
            max number of threads for parts threads are available in this pipeline.
        

        Returns
        -------
        None.

        '''
        super().__init__(reference, fastq, minion, threads)   
        self.reference = SCRIPT_PATH + "refs/HAV.fasta"

    def mapping(self):
        super().mapping()
        bam_path = "BAM/"
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file and "bai" not in bam_file:
                subprocess.call(SPLIT % dict(bam=bam_path + bam_file), shell=True)#split bams by reference
                os.remove(bam_path + bam_file)