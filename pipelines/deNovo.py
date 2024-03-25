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
# RUN_SPADES = "spades --12 %(r1)s -o %(output_path)s --rna"

RUN_BLAST = "blastn -db %(db)s -query %(query)s -out %(output_file)s -outfmt \"6  qseqid sseqid stitle qcovs qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore \""


class de_novo(general_pipe):
      
    def __init__(self, reference, fastq, minion, threads, virus_name):
        super().__init__(reference, fastq, minion, threads) 
        self.virus_name = virus_name if not virus_name==True else ''
        if self.minion:
            raise ValueError("de-novo pipeline only works with illumina reads.")
        create_dirs(["BAM/fastq_based", "BAM/contig_based", "VCF/fastq_based", "VCF/contig_based"]) #temp comment
        
    #run spades on paired end fastq files. the list must contain sample, r1 r2 file names 
    def run_spades(self):
        for sample, r1 in self.sample_fq_dict.items():
            r2 = r1.replace("R1","R2")
            create_dirs(["spades/spades_results/"+sample])#temp comment
            subprocess.call(RUN_SPADES % dict( r1=self.fastq + r1, r2=self.fastq + r2, output_path="spades/spades_results/"+ sample + "/"), shell=True)

    #run blastn on spades contigs - trascript.fasta  
    def run_blast(self):
        create_dirs(["blast"])#temp comment
        for spades_result in os.listdir("spades/spades_results/"):
            query="spades/spades_results/" + spades_result + "/transcripts.fasta"
            blast_result = "blast/" + spades_result + ".txt"
            if os.path.isfile(query):
                subprocess.call(RUN_BLAST % dict(db=self.reference, query=query, output_file="blast/" + spades_result + ".txt"), shell=True)
            else:
                with open(blast_result, 'w'):
                    pass
            # csv format
            dataframe = pd.read_csv(blast_result, delimiter="\t",
                                    names=["query_seq_id","subject_seq_id", "subject_title", 
                                           "query_coverage", "query_length" , "subject_length",
                                           "identity", "alignment_length", "mismatches", "gap_open",
                                           "query_start", "query_end", "subject_start", "subject_end",
                                           "E-value", "bitscore"], index_col=False)
            dataframe.to_csv(blast_result.replace("txt", "csv"), encoding='utf-8', index=False)
            # os.remove(blast_result)
    #load blast output and choose the reference of the longest highest score contig. 
    def choose_reference_filter_contigs(self):
        create_dirs(["fasta/","fasta/selected_contigs","fasta/all_contigs","fasta/selected_references","BAM/fastq_based/","BAM/contig_based/"])
        sample_ref = pd.DataFrame(columns=["sample", "reference_id", "reference_name","blast_coverage(contig_to_ref)", "blast_identity", "blast_score"]) #df of the selected reference for each sample
        for sample, r1 in self.sample_fq_dict.items():
           df = pd.read_csv('blast/'+sample+'.csv')
           if len(df) < 0:
               continue
    
           #filter the highest score of each contig
           df = df[~df['subject_title'].str.contains('retrovirus|endogenous|phage', case=False)]
           prefered_virus = df[df['subject_title'].str.contains(self.virus_name, case=False)]
           prefered_virus_ref =  prefered_virus[prefered_virus['subject_title'].str.contains('complete', case=False)]
           if len(prefered_virus_ref) > 0:
               df = prefered_virus_ref
               
           best_refs = df.groupby("query_seq_id")['bitscore'].max()
           df = df[df['bitscore'].isin(best_refs)]
           df.to_csv('blast/filtered_'+sample+".csv", index=False)
           
           #create multi-fasta contains only the contigs blast found
           contigs_list = df['query_seq_id'].tolist()
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
           longest_contig = df.loc[df['query_length'].idxmax()]
           reference_id = longest_contig['subject_seq_id']
           reference_id = reference_id.split("|")[1].split("|")[0] if "|" in reference_id else reference_id
           reference_name = longest_contig['subject_title']
           blast_score = longest_contig['bitscore']
           ident = longest_contig['identity']
           coverage = longest_contig['query_coverage']
           sample_ref.loc[len(sample_ref)] = [str(sample),reference_id, reference_name, coverage, ident, blast_score]
        #sample_ref.to_csv("/home/hagar/sample_ref.csv", index = False)
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
        self.run_spades()
        self.run_blast()
        self.sample_ref = self.choose_reference_filter_contigs()
        self.import_references()
        
        filter_out_code = 4
        
        #align selected contig 
        for sample, r1 in self.sample_fq_dict.items():
            r2 = r1.replace("R1","R2")
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
    def cns(self, bam_path, cns_path, cns_x_path, min_depth_call, min_freq_thresh):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        create_dirs([cns_path+contig_dir, cns_path+fastq_dir, cns_x_path+contig_dir,cns_x_path+fastq_dir])
        super().cns(bam_path+contig_dir, cns_path+contig_dir, cns_x_path+contig_dir, min_depth_call, min_freq_thresh)
        super().cns(bam_path+fastq_dir,cns_path+fastq_dir, cns_x_path+fastq_dir, min_depth_call, min_freq_thresh)
        return
    
    def depth(self, bam_path, depth_path):
        contig_dir = "contig_based/"
        fastq_dir = "fastq_based/"
        create_dirs([depth_path+contig_dir, depth_path+fastq_dir])
        super().depth(bam_path+contig_dir, depth_path+contig_dir)
        super().depth(bam_path+fastq_dir, depth_path+fastq_dir)
    #override
    def mafft(self, not_aligned, aligned):
        return
    
    #override
    #write report twice (contigs based and fastq based)
    def qc_report(self, bam_path, depth_path, output_report):
        for report_type in ['contig_based', 'fastq_based']:  
            super().qc_report(bam_path+report_type+'/', depth_path+report_type+'/', output_report + "_" + report_type)
        
            #add redference columns 
            report = pd.read_csv(output_report+"_"+ report_type + ".csv", dtype={'sample': str})
            
            report = self.sample_ref.merge(report, how='left', on='sample')
            report.to_csv(output_report + "_" + report_type + '.csv', index = False)    
        
        
        