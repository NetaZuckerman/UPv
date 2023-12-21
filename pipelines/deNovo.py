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
from pipelines.generalPipeline import general_pipe, INDEX, FILTER_BAM, SORT,CHIMER


#alignment conmmands
BWM_MEM_CONTIGS = "bwa mem -v1 -t  %(threads)s  %(reference)s %(sample_fasta)s | samtools view -@ 32 -b - > %(output_path)s%(out_file)s.bam"
BWM_MEM_FASTQ = "bwa mem -v1 -t  %(threads)s  %(reference)s %(r1)s %(r2)s | samtools view -@ 32 -b - > %(output_path)s%(out_file)s.bam"

RUN_SPADES = "spades -1 %(r1)s -2 %(r2)s -o %(output_path)s --rna"
#RUN_SPADES = "spades --12 %(r1)s -o %(output_path)s --rna"

RUN_BLAST = "blastn -db %(db)s -query %(query)s -out %(output_file)s -outfmt \"6 qseqid sseqid stitle qlen qcovs score bitscore \""


class de_novo(general_pipe):
      
    def __init__(self, reference, fastq, threads, virus_name):
        super().__init__(reference, fastq, threads) 
        self.virus_name = virus_name if not virus_name==True else ''
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
        sample_ref = pd.DataFrame(columns=["sample", "reference_id", "reference_name"]) #df of the selected reference for each sample
        for sample, r1, r2 in self.r1r2_list:   
           if os.stat('blast/'+sample+'.txt').st_size == 0:
               continue
           df = pd.read_csv('blast/'+sample+'.txt', sep="\t", header=None)
           df.columns = ["contig_seq-id", "reference_seq-id", "reference_title", "contig_length" ,"coverage(contig_to_ref)", "raw_score", "bit_score"] 
           #filter the highest score of each contig
           df = df[~df['reference_title'].str.contains('retrovirus|endogenous|phage', case=False)]
           prefered_virus = df[df['reference_title'].str.contains(self.virus_name, case=False)]
           prefered_virus_ref =  prefered_virus[prefered_virus['reference_title'].str.contains('complete', case=False)]
           if len(prefered_virus_ref) > 0:
               df = prefered_virus_ref
               
           best_refs = df.groupby("contig_seq-id")['raw_score'].max()
           df = df[df['raw_score'].isin(best_refs)]
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
           reference_id = df.loc[df['contig_length'].idxmax()]['reference_seq-id']
           reference_id = reference_id.split("|")[1].split("|")[0] if "|" in reference_id else reference_id
           reference_name = df.loc[df['contig_length'].idxmax()]['reference_title']
           sample_ref.loc[len(sample_ref)] = [sample,reference_id, reference_name]
        return sample_ref
    
    #export reference sequence fasta 
    def import_references(self):
        for i, row in self.sample_ref.iterrows():
            ref_id = row["reference_id"]
            sample = row["sample"]
            for record in SeqIO.parse(self.reference, "fasta"):
                if ref_id in record.description:
                    ref_file = open("fasta/selected_references/"+sample+".fasta", 'a')
                    SeqIO.write(record,ref_file,"fasta")
                    ref_file.close()
                    continue
    #override
    #map each sample to its reference
    def mapping(self):
        # self.run_spades()
        # self.run_blast()
        self.sample_ref = self.choose_reference_filter_contigs()
        self.import_references()
        
        filter_out_code = 4
        
        #align selected contig 
        for sample, r1, r2 in self.r1r2_list:
            fasta = sample + ".fasta"
            if os.path.exists("fasta/selected_references/" + fasta):
                subprocess.call(INDEX % dict(reference="fasta/selected_references/" + fasta), shell=True)
                subprocess.call(BWM_MEM_CONTIGS % dict(reference="fasta/selected_references/" + fasta, sample_fasta="fasta/selected_contigs/"+fasta, out_file=sample, output_path="BAM/contig_based/",threads = self.threads), shell=True) #map to reference
                subprocess.call(BWM_MEM_FASTQ % dict(reference="fasta/selected_references/" + fasta ,r1=self.fastq+r1, r2=self.fastq+r2, out_file=sample, output_path="BAM/fastq_based/",threads = self.threads), shell=True) #map to reference
                subprocess.call(FILTER_BAM % dict(sample=sample,filter_out_code = filter_out_code, output_path="BAM/fastq_based/",threads = self.threads), shell=True) #keep mapped reads
                subprocess.call(FILTER_BAM % dict(sample=sample,filter_out_code = filter_out_code, output_path="BAM/contig_based/",threads = self.threads), shell=True) #keep mapped reads
                subprocess.call(SORT % dict(sample=sample, output_path="BAM/fastq_based/",threads = self.threads), shell=True)         
                subprocess.call(SORT % dict(sample=sample, output_path="BAM/contig_based/",threads = self.threads), shell=True)   
                subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/fastq_based/"), shell=True)
                subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/contig_based/"), shell=True)
    
    #override
    #run general pipeline function twice (contigs based and fastq based)
    def cns_depth(self, bam_path, depth_path, cns_path, cns5_path, thresh):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        create_dirs([depth_path+contig_dir, depth_path+fastq_dir, cns_path+contig_dir, cns_path+fastq_dir, cns5_path+contig_dir,cns5_path+fastq_dir])
        super().cns_depth(bam_path+contig_dir, depth_path+contig_dir, cns_path+contig_dir, cns5_path+contig_dir, thresh)
        super().cns_depth(bam_path+fastq_dir, depth_path+fastq_dir, cns_path+fastq_dir, cns5_path+fastq_dir, thresh)
        return
    
    #override
    def mafft(self, not_aligned, aligned):
        return
    
    #override
    #write report twice (contigs based and fastq based)
    def qc_report(self, bam_path, depth_path, output_report, vcf=0):
        for report_type in ['contig_based', 'fastq_based']:  
            super().qc_report(bam_path+report_type+'/', depth_path+report_type+'/', output_report + "_" + report_type)
        
            #add redference columns 
            report = pd.read_csv(output_report+"_"+ report_type + ".csv")
            report = self.sample_ref.merge(report, how='left', on='sample')
            report.to_csv(output_report + "_" + report_type + '.csv', index = False)    
        
        
        