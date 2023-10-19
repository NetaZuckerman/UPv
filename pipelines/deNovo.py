#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:29:54 2022

@author: hagar

this script gets fastq path of paired end read with no reference.
it runs rna SPAdes and saves the output for each sample in spades/spades_results.
in rnaSPAdes the contigs + scalffolds are in "transcript.fasta".

"""
from utils.utils import create_dirs
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import shutil
from pipelines.generalPipeline import general_pipe, INDEX, FILTER_BAM, SORT


#alignment conmmands
BWM_MEM_CONTIGS = "bwa mem -v1 -t 32 %(reference)s %(sample_fasta)s | samtools view -@ 32 -b - > %(output_path)s%(out_file)s.bam"
BWM_MEM_FASTQ = "bwa mem -v1 -t 32 %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b - > %(output_path)s%(out_file)s.bam"

RUN_SPADES = "spades -1 %(r1)s -2 %(r2)s -o %(output_path)s --rna"
#RUN_SPADES = "spades --12 %(r1)s -o %(output_path)s --rna"

RUN_BLAST = "blastn -db %(db)s -query %(query)s -out %(output_file)s -outfmt \"6 qseqid sseqid stitle qlen qcovs score bitscore \""


class de_novo(general_pipe):
      
    def __init__(self, reference, fastq):
        super().__init__(reference, fastq) 
        create_dirs(["BAM/fastq_based", "BAM/contig_based", "VCF/fastq_based", "VCF/contig_based"]) #temp comment
        
    #run spades on paired end fastq files. the list must contain sample, r1 r2 file names 
    def run_spades(self):
        for sample, r1, r2 in self.r1r2_list:
            create_dirs(["spades/spades_results/"+sample])#temp comment
            subprocess.call(RUN_SPADES % dict( r1=self.fastq + r1, r2=self.fastq + r2, output_path="spades/spades_results/"+ sample + "/"), shell=True)


    #run blastn on spades contigs - trascript.fasta  
    def run_blast(self):
        create_dirs(["blast"])#temp comment
        for spades_result in os.listdir("spades/spades_results/"):
            subprocess.call(RUN_BLAST % dict(db=self.reference, query="spades/spades_results/" + spades_result + "/transcripts.fasta", output_file="blast/" + spades_result + ".txt"), shell=True)
            
            # csv format
            dataframe = pd.read_csv("blast/" + spades_result + ".txt",delimiter="\t", names=["query_sequence", "query_seq_id", "subject_title" ,"query_length" ,"query_coverage", "raw_score", "bit_score"])
            dataframe.to_csv("blast/" + spades_result + ".csv", encoding='utf-8', index=False)

    #load blast output and choose the reference of the longest highest score contig. 
    def choose_reference_filter_contigs(self):
        create_dirs(["fasta/","fasta/selected_contigs","fasta/all_contigs","fasta/selected_references","BAM/fastq_based/","BAM/contig_based/"])
        sample_ref = {} #dict of the selected reference for each sample
        for sample, r1, r2 in self.r1r2_list:   
           if os.stat('blast/'+sample+'.txt').st_size == 0:
               continue
           df = pd.read_csv('blast/'+sample+'.txt', sep="\t", header=None)
           df.columns = ["contig_seq-id", "reference_seq-id", "reference_title", "contig_length" ,"coverage(contig_to_ref)", "raw_score", "bit_score"] 
           #filter the highest score of each contig
           approved_fields = df.groupby("contig_seq-id")['raw_score'].max()
           df = df[df['raw_score'].isin(approved_fields)]
           df.to_csv('blast/filtered_'+sample+".csv", index=False)
           
           #create multi-fasta contains only the contigs blast found
           contigs_list = df['contig_seq-id'].tolist()
           #write fasta file with the selected contigs
           filtered_file = open("fasta/selected_contigs/"+sample+".fasta", 'a')
           #move spades contigs 
           shutil.copyfile("spades/spades_results/"+sample+"/transcripts.fasta","fasta/all_contigs/" +sample+".fasta")
           
           for record in SeqIO.parse("fasta/all_contigs/"+sample+".fasta", "fasta"):
               if record.description in contigs_list:
                   filtered_file.write(record.format("fasta"))
           filtered_file.close()
       #shutil.rmtree("spades/spades_results/")
            #select the reference of the longest contig
           reference = df.loc[df['contig_length'].idxmax()]['reference_seq-id'].split("|")[1].split("|")[0]
           sample_ref.update({sample:reference})
        return sample_ref
    
    #export reference sequence fasta 
    def import_references(self, sample_ref):
        for sample, ref in sample_ref.items():
            for record in SeqIO.parse(self.reference, "fasta"):
                if ref in record.description:
                    ref_file = open("fasta/selected_references/"+sample+".fasta", 'a')
                    SeqIO.write(record,ref_file,"fasta")
                    ref_file.close()
                    continue
    #override
    #map each sample to its reference
    def bam(self,sample_r1_r2):
        filter_out_code = 4
        #align selected contig 
        sample = sample_r1_r2[0]
        r1= sample_r1_r2[1]
        r2 = sample_r1_r2[2]
        fasta = sample + ".fasta"
        if os.path.exists("fasta/selected_references/" + fasta):
            subprocess.call(INDEX % dict(reference="fasta/selected_references/" + fasta), shell=True)
            subprocess.call(BWM_MEM_CONTIGS % dict(reference="fasta/selected_references/" + fasta, sample_fasta="fasta/selected_contigs/"+fasta, out_file=sample, output_path="BAM/contig_based/"), shell=True) #map to reference
            subprocess.call(BWM_MEM_FASTQ % dict(reference="fasta/selected_references/" + fasta ,r1=self.fastq+r1, r2=self.fastq+r2, out_file=sample, output_path="BAM/fastq_based/"), shell=True) #map to reference
            subprocess.call(FILTER_BAM % dict(sample=sample,filter_out_code = filter_out_code, output_path="BAM/fastq_based/"), shell=True) #keep mapped reads
            subprocess.call(FILTER_BAM % dict(sample=sample,filter_out_code = filter_out_code, output_path="BAM/contig_based/"), shell=True) #keep mapped reads
            subprocess.call(SORT % dict(sample=sample, output_path="BAM/fastq_based/"), shell=True)         
            subprocess.call(SORT % dict(sample=sample, output_path="BAM/contig_based/"), shell=True)         
    
    #override
    #run general pipeline function twice (contigs based and fastq based)
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        create_dirs([depth_path+contig_dir, depth_path+fastq_dir, cns_path+contig_dir, cns_path+fastq_dir, cns5_path+contig_dir,cns5_path+fastq_dir])
        super().cns_depth(bam_path+contig_dir, depth_path+contig_dir, cns_path+contig_dir, cns5_path+contig_dir)
        super().cns_depth(bam_path+fastq_dir, depth_path+fastq_dir, cns_path+fastq_dir, cns5_path+fastq_dir)
    
    #override
    def variant_calling(self, bam_path, vcf_path):
        super().variant_calling(bam_path = "BAM/contig_based/", vcf_path= "VCF/contig_based/")
        super().variant_calling(bam_path = "BAM/fastq_based/", vcf_path= "VCF/fastq_based/")
    
    #override
    #write report twice (contigs based and fastq based)
    def qc_report(self, bam_path, depth_path, output_report, vcf=0):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"        
        super().results_report(bam_path+contig_dir, depth_path+contig_dir, output_report+"_contig_based")
        super().results_report(bam_path+fastq_dir, depth_path+fastq_dir, output_report+"_fastq_based")
        
        

        
        
        
        