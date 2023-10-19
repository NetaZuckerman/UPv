#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 11:16:08 2022

@author: hagar

Additional analysis for polio data.

"""


RUN_SPADES = "spades --12 %(sample)s -o %(output_path)s --rna"
FILTER = "bwa mem -v1 -t 16 %(reference)s %(r1)s %(r2)s | samtools view -@ 16 -b -f 2 > %(output_path)s%(sample)s.bam"
min_FILTER = "bwa mem -v1 -t 16 %(reference)s %(fastq)s | samtools view -@ 16 -b -F 4 > %(output_path)s%(sample)s.bam"
BAM2FQ = "samtools bam2fq -0n %(file)s.bam > %(file)s.fastq"
minion_BAM2FQ = "samtools bam2fq %(file)s.bam > %(file)s.fastq"
#BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #filter out common reads
BWM_MEM_FASTQ = "bwa mem -v1 -t 16 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 16 -b > %(output_path)s%(out_file)s.bam" #filter out common reads
BWM_MEM_CONTIGS = "bwa mem -v1 -t 16  %(reference)s %(sample_fasta)s | samtools view -bq 1 | samtools view -@ 16 -b -F 4- > %(output_path)s%(out_file)s.bam"
#old mappings:###
    #BWM_MEM_FASTQ = "bwa mem -v1 -t 32 -B 10000 %(reference)s %(fastq)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #strict mapping (high penalty)
    #BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(fastq)s | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam" #no filters
    #BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -bq 1 | samtools view -@ 32 -b -F 4- > %(output_path)s%(out_file)s.bam"
#################
from utils import utils
from pipelines.generalPipeline import  FILTER_BAM, SORT, general_pipe
import subprocess
import os
from utils.utils import SPLIT
from utils.summerize_coverage import summerize_coverage

class polio(general_pipe):
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq) 
        utils.create_dirs([self.fastq+"polio_reads","BAM/fastq_based", "BAM/contig_based"]) #temp comment

    #filter the fastq files to contain only mapped reads 
    #in minion the we dont have r1 and r2. all files should be merged by barcodes before running this code
    #TODO: implement minion

    def filter_not_polio(self):
        '''
        filter outs reads that are not mapped to polio.
        
        Parameters
        ----------

        Returns
        -------
        None.

        '''
        for sample, r1, r2 in self.r1r2_list:
            subprocess.call(FILTER % dict( reference=self.reference, r1=self.fastq + r1, r2=self.fastq + r2, output_path=self.fastq+"polio_reads/", sample=sample), shell=True)
            subprocess.call(BAM2FQ % dict(file=self.fastq+"polio_reads/" + sample), shell=True)
            os.remove(self.fastq+"polio_reads/" + sample + ".bam")
    
    def run_spades(self):
        for sample, r1, r2 in self.r1r2_list:
            utils.create_dirs(["spades/spades_results/"+sample])
            subprocess.call(RUN_SPADES % dict( sample= self.fastq+sample+".fastq", output_path="spades/spades_results/"+ sample + "/"), shell=True)
              
    #override
    #map each sample to its reference
    def bam(self,treads):
        for sample, r1, r2 in self.r1r2_list:
            subprocess.call(BWM_MEM_CONTIGS % dict(reference=self.reference, sample_fasta="spades/spades_results/"+ sample + "/"+"transcripts.fasta", out_file=sample, output_path="BAM/contig_based/"), shell=True) #map to reference
            subprocess.call(BWM_MEM_FASTQ % dict(reference=self.reference ,fastq=self.fastq+sample+".fastq", out_file=sample, output_path="BAM/fastq_based/"), shell=True) #map to reference
            subprocess.call(SPLIT % dict(bam="BAM/fastq_based/"+sample+".bam"), shell=True)
            os.remove("BAM/fastq_based/"+sample+".bam")
            subprocess.call(SPLIT % dict(bam="BAM/contig_based/"+sample+".bam"), shell=True)
            os.remove("BAM/contig_based/"+sample+".bam")
            self.map_bam()

    def map_bam(self):
        filter_out_code = 4
        for bam in os.listdir("BAM/fastq_based/"):
            sample_ref = bam.split(".bam")[0]
            subprocess.call(FILTER_BAM % dict(sample=sample_ref, filter_out_code = filter_out_code, output_path="BAM/fastq_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample_ref, output_path="BAM/fastq_based/"), shell=True)         
        for bam in os.listdir("BAM/contig_based/"):
            sample_ref = bam.split(".bam")[0]
            subprocess.call(FILTER_BAM % dict(sample=sample_ref, filter_out_code = filter_out_code, output_path="BAM/contig_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample_ref, output_path="BAM/contig_based/"), shell=True)         
        
        
    #override
    #run general pipeline function twice (contigs based and fastq based)
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path, cnsThresh):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        utils.create_dirs([depth_path+contig_dir, depth_path+fastq_dir, cns_path+contig_dir, cns_path+fastq_dir, cns5_path+contig_dir,cns5_path+fastq_dir])
        super().cns_depth(bam_path+contig_dir, depth_path+contig_dir, cns_path+contig_dir, cns5_path+contig_dir, cnsThresh)
        super().cns_depth(bam_path+fastq_dir, depth_path+fastq_dir, cns_path+fastq_dir, cns5_path+fastq_dir, cnsThresh)
       
    #override TODO - tests
    def variant_calling(self, bam_path, vcf_path):
        utils.create_dirs(["VCF/fastq_based", "VCF/contig_based"])
        super().variant_calling(bam_path = "BAM/contig_based/", vcf_path= "VCF/contig_based/")
        super().variant_calling(bam_path = "BAM/fastq_based/", vcf_path= "VCF/fastq_based/")
        
    #override
    #write report twice (contigs based and fastq based)
    def qc_report(self, bam_path, depth_path, output_report, vcf=0):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"        
        super().results_report(bam_path+contig_dir, depth_path+contig_dir, output_report+"_contig_based") # TODO - fix de novo report
        super().results_report(bam_path+fastq_dir, depth_path+fastq_dir, output_report+"_fastq_based")
        
        summerize_coverage(output_report+"_fastq_based.csv")
        summerize_coverage(output_report+"_contig_based.csv")

        
        
        
        
        
        
        
        
        
        
        
    