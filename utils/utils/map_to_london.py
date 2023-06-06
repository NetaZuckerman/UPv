#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 12:18:52 2022

@author: hagar

Filter sabin 1 and sabin3 reads and map to london

"""
filter_path = "/mnt/project1/projects/POLIO/VDPV2_2022/tests_270722/london_wg/filtered_fastq/"
MAP1 = "bwa mem -v1 -t 32 -B 10000 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b -f 4- > "+ filter_path +"%(sample)s.s1.unmapped.strict.bam"
MAP3 = "bwa mem -v1 -t 32 -B 10000 %(reference)s "+ filter_path +"%(fastq)s| samtools view -@ 32 -b -f 4- >"+ filter_path +"%(sample)s.s3.unmapped.strict.bam"
bam2fq = "samtools bam2fq "+ filter_path +"%(sample)s.%(reference)s.unmapped.strict.bam > "+ filter_path +"%(sample)s.%(reference)s.unmapped.strict.fastq"
bam2fq2 = "samtools bam2fq "+ filter_path +"%(sample)s.%(reference)s.unmapped.strict.bam > "+ filter_path +"%(sample)s.fastq"
fastq_path= "/mnt/project1/projects/POLIO/VDPV2_2022/fastq/"
filter_path = "/mnt/project1/projects/POLIO/VDPV2_2022/tests_270722/london_wg/filtered_fastq/"
bam_path = "/mnt/project1/projects/POLIO/VDPV2_2022/tests_270722/london_wg/BAM/"
s1 = "/home/hagar/refs/polio/Sabin1_AY184219.1.fasta"
s3 = "/home/hagar/refs/polio/Sabin3.fasta"
london = "/mnt/project1/projects/POLIO/VDPV2_2022/london/WGS/our_vdpv2_mapped_to_london_whole_genome/london_wg.fasta"

from utils import get_r1r2_list
import subprocess
r1_r2_list = get_r1r2_list(fastq_path)
for r1_r2 in r1_r2_list:
    sample = r1_r2[0]
    r1 = fastq_path + r1_r2[1]
    r2 = fastq_path + r1_r2[2]
        
    subprocess.call(MAP1 % dict(reference=s1, r1=r1, r2=r2 , sample=sample), shell=True)
    subprocess.call(bam2fq % dict(reference="s1", sample=sample), shell=True)
    subprocess.call(MAP3 % dict(reference=s3, fastq=sample+".s1.unmapped.strict.fastq" , sample=sample), shell=True)
    subprocess.call(bam2fq2 % dict(reference="s3", sample=sample), shell=True)


















