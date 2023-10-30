#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 08:05:15 2023

@author: hagar
"""

#regions
ul30_reg = (62807, 66512) 	
ul42_reg = (93113, 94577)
ul23_reg = (46675, 47803)


import pandas as pd
from pipelines.generalPipeline import general_pipe
import os 
from utils.utils import get_sequences, write_sub_fasta, create_dirs
from mutations import signatures
from utils.format_xl import save_format_xl


SCRIPT_PATH = os.path.dirname(__file__) +'/'

class hsv(general_pipe):
    def __init__(self, reference, fastq, threads):
        super().__init__(reference, fastq, threads)    
 

#cut by gene
    def cut_genes(self, aln_path):
        fasta = get_sequences(aln_path + "all_aligned.fasta")
        write_sub_fasta(fasta, aln_path, ul30_reg, "reg_UL30")
        write_sub_fasta(fasta, aln_path, ul42_reg, "reg_UL42")
        write_sub_fasta(fasta, aln_path, ul23_reg, "reg_UL23", '-')
        
#mutations
    def mutations(self):
        create_dirs(["reports"]) 
        regions_file_UL23 = SCRIPT_PATH + "refs/herpes1_regions_UL23.csv"
        regions_file_UL30 = SCRIPT_PATH + "refs/herpes1_regions_UL30.csv"
        regions_file_UL42 = SCRIPT_PATH + "refs/herpes1_regions_UL42.csv"
        signatures.run("alignment/reg_UL23.fasta", regions_file_UL23, "reports/mutations_UL23.xlsx")
        signatures.run("alignment/reg_UL30.fasta", regions_file_UL30, "reports/mutations_UL30.xlsx")
        signatures.run("alignment/reg_UL42.fasta", regions_file_UL42, "reports/mutations_UL42.xlsx")


    def resist_poly_reports(self):
        '''
        use the mutations table to generate resistant and polymorphism reports

        Returns
        -------
        None.

        '''
        
        ####resistant mutations - SNPs and deletions
        resist = pd.read_excel(SCRIPT_PATH + "refs/hsv_resist_mut.xlsx", "SNP&del")
        mut_tbl_UL23 = pd.read_excel("reports/mutations_UL23.xlsx")
        mut_tbl_UL23["nt_position_on_genome"] = mut_tbl_UL23["nt_position_on_genome"] + ul23_reg[0]-1
        mut_tbl_UL30 = pd.read_excel("reports/mutations_UL30.xlsx")
        mut_tbl_UL30["nt_position_on_genome"] = mut_tbl_UL30["nt_position_on_genome"] + ul30_reg[0]-1
        mut_tbl_UL42 = pd.read_excel("reports/mutations_UL42.xlsx")
        mut_tbl_UL42["nt_position_on_genome"] = mut_tbl_UL42["nt_position_on_genome"] + ul42_reg[0]-1
        
        mut_tbl = pd.concat([mut_tbl_UL23,mut_tbl_UL30,mut_tbl_UL42])
        merged = mut_tbl.merge(resist, how='left', on=["gene_name","nt_position_on_gene" ,"X14112.1_NT"])
        seq_num = int((len(merged.columns) - len(resist.columns)) /2 ) -3
        save_format_xl(merged, seq_num, "reports/mutations_snp&del.xlsx")
        
        #######resistant mutations - insertions
        res_inser = pd.read_excel(SCRIPT_PATH + "refs/hsv_resist_mut.xlsx", "insertions")
        inser = pd.read_csv("alignment/all_aligned.fasta.insertions.csv")
        
        #reshape insertions df to match resist_insertions db and extrant insertions positions
        inser = inser.set_index("strain").transpose().reset_index().rename(columns={"index" : "nt_position_on_genome"}) 
        inser["nt_position_on_genome"] = inser["nt_position_on_genome"].str.split(" ").str[-1].astype({'nt_position_on_genome': 'int32'})
        merged = inser.merge(res_inser, how='left', on=["nt_position_on_genome"])
        merged.to_csv("reports/mutations_insertions.csv", index=False)
        
        ####polymorphism nucleotide
        poly = pd.read_csv(SCRIPT_PATH + "refs/hsv_polymorphism_nucleotides.csv")
        merged = mut_tbl.merge(poly, how='left', on=["gene_name","nt_position_on_gene" ,"X14112.1_NT"]).dropna()
        save_format_xl(merged, seq_num, "reports/polymorphism.xlsx")

        #remove script temp files:
        os.remove("reports/mutations_UL23.xlsx")
        os.remove("reports/mutations_UL30.xlsx")
        os.remove("reports/mutations_UL42.xlsx")            
    
        
#override report
    def qc_report(self, bam_path, depth_path, output_report):
        self.cut_genes("alignment/")
        self.mutations()
        self.resist_poly_reports()
        
        general_qc = pd.DataFrame(columns=['sample', 'mapped%','mapped_reads','total_reads','cov_bases',"chimeric_read_count"])
        ul23_qc = pd.DataFrame(columns=['sample','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth'])
        ul30_qc = pd.DataFrame(columns=['sample','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth'])
        ul42_qc = pd.DataFrame(columns=['sample','coverage%','coverage_CNS_5%', 'mean_depth','max_depth','min_depth'])
        
        for bam_file in os.listdir(bam_path):
                if "sorted" in bam_file and "bai" not in bam_file:
                    #general qc
                    sample = bam_file.split(".mapped")[0] + bam_file.split(".sorted")[1].split(".bam")[0]
                    total_reads, mapped_reads, mapped_percentage, cov_bases, chimer_count = super().general_qc(bam_path+bam_file)
            
                    #depth 
                    max_depth23, min_depth23,mean_depth23, cover23, cns5_cover23 = super().depth(depth_path+sample+".txt", start = ul23_reg[0], end= ul23_reg[1])
                    max_depth30, min_depth30, mean_depth30, cover30, cns5_cover30 = super().depth(depth_path+sample+".txt", start = ul30_reg[0], end= ul30_reg[1])
                    max_depth42, min_depth42,mean_depth42, cover42, cns5_cover42 = super().depth(depth_path+sample+".txt", start = ul42_reg[0], end= ul42_reg[1])
                    
                    #add rows
                    general_qc.loc[len(general_qc.index)] = [sample, mapped_percentage, mapped_reads, total_reads, cov_bases, chimer_count]
                    ul23_qc.loc[len(ul23_qc.index)] = [sample, cover23, cns5_cover23, mean_depth23, max_depth23, min_depth23]
                    ul30_qc.loc[len(ul30_qc.index)] = [sample, cover30, cns5_cover30, mean_depth30, max_depth23, min_depth30]
                    ul42_qc.loc[len(ul42_qc.index)] = [sample, cover42, cns5_cover42, mean_depth42, max_depth42, min_depth42]
                    
        with pd.ExcelWriter('QC/QC_report.xlsx') as writer:  
            general_qc.to_excel(writer, sheet_name='general')
            ul23_qc.to_excel(writer, sheet_name='UL23')
            ul30_qc.to_excel(writer, sheet_name='UL30')
            ul42_qc.to_excel(writer, sheet_name='UL42')
                        
         