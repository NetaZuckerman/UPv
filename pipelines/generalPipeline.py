#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:37:02 2022

@author: hagar

General pipeline is the base class of the upv pipeline.
"""

import subprocess
import os
import pysam
from statistics import mean
import csv
from utils import utils
import shutil

MAIN_SCRIPT_DIR = os.path.dirname(__file__)+'/../'


#alignment conmmands
INDEX = "bwa index %(reference)s"
BWM_MEM = "bwa mem -v1 -t %(threads)s %(reference)s %(r1)s %(r2)s | samtools view -@ %(threads)s -b - > %(output_path)s%(sample)s.bam"
MINIMAP = "minimap2 -ax map-ont %(ref)s %(fastq_dir)s/*.fastq* > %(output)s.bam"
FILTER_BAM = "samtools view -@ %(threads)s -b -F %(filter_out_code)s %(output_path)s%(sample)s.bam > %(output_path)s%(sample)s.mapped.bam"
SORT = "samtools sort -@ %(threads)s %(output_path)s%(sample)s.mapped.bam -o %(output_path)s%(sample)s.mapped.sorted.bam"
CHIMER = "samtools view %(output_path)s%(sample)s.bam |  grep 'SA:' > %(output_path)s%(sample)s.chimeric_reads.txt"
DEPTH = "samtools depth -a %(bam_path)s%(bam_file)s > %(depth_path)s%(sample)s.txt"
CNS = "samtools mpileup -A %(bam_path)s%(bam_file)s | ivar consensus -t %(min_freq_thresh)s -m %(min_depth_call)s -p %(cns_path)s%(sample)s.fa"

#MAFFT commands
ALL_NOT_ALIGNED =   "cat %(dir)s > alignment/all_not_aligned.fasta"

#report commands
BREADTH_CNS5 = "$(cut -f3 QC/depth/%(sample)s.txt | awk '$1>5{c++} END{print c+0}')"
#directories for output

            
class general_pipe():
    
    
    def __init__(self, reference, fastq, minion,threads):
        self.reference = reference
        self.fastq = fastq
        self.minion = minion
        if not minion:
            self.sample_fq_dict = utils.get_sample_fq_dict(self.fastq)
        else:
            self.sample_fq_dict = utils.get_barcodes(minion)
        self.threads = threads


    def mapping(self):
        '''
        generate bam file from paired-end fastq (R1,R2).
        generate filtered file with mapped reads.
        generate sorted file
        
        Parameters
        ----------
        sample_r1_r2 : list of lists
            [sample, r1 fastq path, r2 fastq path].

        '''
        subprocess.call(INDEX % dict(reference=self.reference), shell=True)

        filter_out_code = 2052
        
        for sample, fq in self.sample_fq_dict.items():
            if self.minion:
                barcode = fq
                subprocess.call(MINIMAP % dict(ref=self.reference, fastq_dir=self.fastq+barcode,\
                                               output="BAM/"+sample), shell=True)
            else:
                r1 = fq
                r2 = r1.replace("R1","R2")
                subprocess.call(BWM_MEM % dict(threads = self.threads, reference=self.reference,r1=self.fastq + r1, r2=self.fastq + r2, sample=sample, output_path="BAM/"), shell=True) #map to reference
            
            subprocess.call(FILTER_BAM % dict(threads = self.threads, sample=sample, filter_out_code = filter_out_code, output_path="BAM/"), shell=True) #keep reads according to the filter code: https://broadinstitute.github.io/picard/explain-flags.html
            subprocess.call(SORT % dict(threads = self.threads, sample=sample, output_path="BAM/"), shell=True)         
            #chimeric reads: find reads that where mapped to more then one region in the genome
            subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/"), shell=True)
        
    def variant_calling(self, bam_path, vcf_path):
        '''

        Parameters
        ----------
        bam_path : str, optional
            The default is "BAM/".
        vcf_path : str, optional
            The default is "VCF/".

        '''
        
        for bam in os.listdir(bam_path):
            if "sorted" in bam:
                sample = bam.split(".bam")[0].replace(".mapped.sorted", "")
                vcf_file = vcf_path + sample + ".csv"

                bamfile = pysam.AlignmentFile(bam_path + bam, "rb")
                reference = pysam.FastaFile(self.reference)
                with open(vcf_file, "w", newline='') as output_file:
                    fieldnames = ['position', 'ref', 'A', 'G', 'C', 'T', 'N', 'D', 'I', 'depth', '%A', '%G', '%C', '%T', '%N', '%D', '%I']
                    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
                    writer.writeheader()

                    for pileupcolumn in bamfile.pileup(min_base_quality =0):
                        position = pileupcolumn.reference_pos
                        reference_base = reference.fetch(region=bamfile.get_reference_name(0), start=position, end=position+1)
                        nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'D': 0, 'I': 0}  # Initialize counts
                        depth = pileupcolumn.n
                        deletions = 0
                        insertions = 0 

                        for pileupread in pileupcolumn.pileups:
                            if pileupread.is_del:
                                deletions += 1
                            elif pileupread.indel > 0:
                                insertions += 1
                            elif not pileupread.is_refskip:
                                base = pileupread.alignment.query_sequence[pileupread.query_position]
                                if base in nucleotide_counts:
                                    nucleotide_counts[base] += 1
                                else:
                                    nucleotide_counts['N'] += 1  # Count ambiguous bases as 'N'

                        row = {
                            'position': position+1,
                            'ref': reference_base,
                            'A': nucleotide_counts['A'],
                            'G': nucleotide_counts['G'],
                            'C': nucleotide_counts['C'],
                            'T': nucleotide_counts['T'],
                            'N': nucleotide_counts['N'],
                            'D': deletions,
                            'I': insertions,
                            'depth': depth,
                            '%A': round((nucleotide_counts['A'] / depth) * 100, 2) if depth > 0 else 0,
                            '%G': round((nucleotide_counts['G'] / depth) * 100, 2) if depth > 0 else 0,
                            '%C': round((nucleotide_counts['C'] / depth) * 100, 2) if depth > 0 else 0,
                            '%T': round((nucleotide_counts['T'] / depth) * 100, 2) if depth > 0 else 0,
                            '%N': round((nucleotide_counts['N'] / depth) * 100, 2) if depth > 0 else 0,
                            '%D': round((deletions / depth) * 100, 2) if depth > 0 else 0,
                            '%I': round((insertions / depth) * 100, 2) if depth > 0 else 0,
                        }
                        writer.writerow(row)

                bamfile.close()
                reference.close()
                
                
    #find mapping depth and consensus sequence 
    def cns(self, bam_path, cns_path, cns_x_path, min_depth_call, min_freq_thresh):
        '''
        Generate consensus sequences and calculate aligning depths from a bam file.
        min_depth_call: ivar consensus -m is  Minimum depth to call consensus. 
        in this pipeline minimum depth 1 is always generated and addional minimum depth is avaialble in the CNS_X directory (default = 5)
        min_freq_thresh: ivar consensus -t    Minimum frequency threshold(0 - 1) to call consensus. defualt = 0.6
        '''
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file and "bai" not in bam_file:
                sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                #consensus
                #CNS - min depth to call = 1
                subprocess.call(CNS % dict(bam_path=bam_path, bam_file=bam_file, cns_path=cns_path,
                                           sample=sample, min_freq_thresh=min_freq_thresh, min_depth_call=1), shell=True) 
                
                #CNS - min depth to call = X (default = 5)
                subprocess.call(CNS % dict(bam_path=bam_path, bam_file=bam_file, 
                                           cns_path=cns_x_path, sample=sample, min_freq_thresh=min_freq_thresh,
                                           min_depth_call=min_depth_call), shell=True)
               
                #remove qual files
                os.remove(cns_path+sample+".qual.txt")
                os.remove(cns_x_path+sample+".qual.txt")
                utils.fix_cns_header(cns_path)
                utils.fix_cns_header(cns_x_path)
    
    def depth(self, bam_path, depth_path):
        for bam_file in os.listdir(bam_path):
            if "sorted" in bam_file and "bai" not in bam_file:
                sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                subprocess.call(DEPTH % dict(bam_path=bam_path, bam_file=bam_file, depth_path=depth_path, sample=sample), shell=True) 
    
    def mafft(self, not_aligned, aligned):
        '''
        multi-fasta align.
        cat all consensus fasta sequences and run MAFFT. the implementation of MAFFT is in utils.

        '''
        subprocess.call(ALL_NOT_ALIGNED % dict(dir="CNS/*"), shell=True)
        utils.mafft(self.reference, not_aligned, aligned)
    
    #write report.csv - mapping analysis
    #de novo flag was created to generate different output for de novo analysis, used in denovo class and polio class.

    def depth_qc(self, file, start=0, end=-1):
        depths = [int(x.split('\t')[2]) for x in open(file).readlines()]
        
        #set
        if end == -1:
            genome_size = len(depths)   
        else:
            depths = depths[start:end]
            genome_size = end-start
        #calculate depths
        depths = [i for i in depths if i != 0]
        if not depths:
            return '','','','',''
        max_depth = max(depths)
        min_depth = min(depths)
        mean_depth = str(round(mean(depths),3))
        coverage = round(len(depths) / genome_size * 100,3)  if genome_size else ''
        
        #cns 5: position with depth > 5
        breadth_cns5 = len([i for i in depths if i > 5])
        cns5_cover = round(breadth_cns5 / genome_size * 100,3)  if genome_size else ''
        
        return max_depth, min_depth,mean_depth, coverage, cns5_cover
    
    def chimer_count(self,file_name):
        #count chimeric reads
        with open(file_name) as f:
            chimer_count = len(f.readlines())
        return chimer_count
    
    def general_qc(self, bam_file, total_reads):
        coverage_stats = pysam.coverage(bam_file).split("\t")
        mapped_reads = int(coverage_stats[11])
        mapped_percentage = round(mapped_reads/total_reads*100,4) if total_reads else ''
        cov_bases =  int(coverage_stats[12])
        chimer_count = self.chimer_count(bam_file.split(".mapped")[0].split(".REF")[0] + ".chimeric_reads.txt")
        
        return mapped_reads, mapped_percentage, cov_bases, chimer_count
        
    def qc_report(self, bam_path, depth_path, output_report):
        f = open(output_report+".csv", 'w')
        writer = csv.writer(f)
        writer.writerow(['sample', 'mapped%','mapped_reads','total_reads','cov_bases','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth', 'chimeric_read_count'])
        
        prev_sample = ''
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    
                    #general qc
                    if not prev_sample == sample.split(".REF")[0]: #when mapping to multi reference dont calculate again the total reads every iteration
                        total_reads = pysam.AlignmentFile(bam_path + bam_file.split(".mapped")[0]+".bam").count(until_eof=True) 
                        prev_sample = sample.split(".REF")[0]
                        
                    mapped_reads, mapped_percentage, cov_bases, chimer_count = self.general_qc(bam_path+bam_file, total_reads)
                    
                    #depth 
                    max_depth, min_depth, mean_depth, cover, cns5_cover = self.depth_qc(depth_path+sample+".txt")
                                        
                    writer.writerow([sample, mapped_percentage, mapped_reads, total_reads, cov_bases, cover, cns5_cover, mean_depth, max_depth, min_depth, chimer_count])
                
        f.close()
        # self.pdf_report(depth_path, output_report+".csv")
    
    def pdf_report(self, depth_path, qc_path):
        rmd_path = os.path.join(MAIN_SCRIPT_DIR, 'utils', 'report.Rmd')
        shutil.copy(rmd_path, os.getcwd() + "/report.Rmd")
        rscript_command = 'Rscript -e "rmarkdown::render(\'{}\', params = list(depth_path = \'{}\', qc_path = \'{}\'))"'.format("report.Rmd", depth_path, qc_path)
        subprocess.run(rscript_command, shell=True, check=True)
        os.remove(os.getcwd() + "/report.Rmd")
        