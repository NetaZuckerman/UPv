#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 09:31:54 2025

@author: hagar

*no minion option*
"""
from pipelines.generalPipeline import general_pipe


class covid(general_pipe):
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