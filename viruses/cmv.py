#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 05:42:03 2024

@author: hagar
"""

import os 
SCRIPT_PATH = os.path.dirname(__file__) +'/'
from utils.utils import create_dirs
from utils import utils
import subprocess
from pipelines.generalPipeline import general_pipe, BWM_MEM, INDEX, FILTER_BAM, SORT,CHIMER
from mutations import signatures
import pandas as pd
from xlsxwriter.utility import xl_col_to_name
import numpy as np

CAT = "cat %(path)s > %(output)s"


class cmv(general_pipe):
    def __init__(self, reference, fastq, minion, threads):
        super().__init__(reference, fastq, minion, threads)    
        self.references_dict = {'ul54': SCRIPT_PATH + "refs/cmv_FJ527563.1_UL54.fasta", 
                          'ul97': SCRIPT_PATH + "refs/cmv_FJ527563.1_UL97.fasta"}
        create_dirs(['BAM/ul54','BAM/ul97','VCF/ul54', 'VCF/ul97'])
        
    def mapping(self):
        for ref_name, ref_path  in self.references_dict.items():
            
            subprocess.call(INDEX % dict(reference=ref_path), shell=True)
            
            for sample, fq in self.sample_fq_dict.items():
                r1 = fq
                r2 = r1.replace("R1","R2")
                
                subprocess.call(BWM_MEM % dict(threads = self.threads, reference=ref_path,
                                               r1=self.fastq + r1, r2=self.fastq + r2, 
                                               sample=sample, output_path="BAM/" + ref_name + '/'), shell=True) #map to reference
                
                subprocess.call(FILTER_BAM % dict(threads = self.threads, sample=sample,
                                                  filter_out_code = 4,
                                                  output_path="BAM/"+ ref_name + '/'), shell=True) 
                
                subprocess.call(SORT % dict(threads = self.threads, sample=sample,
                                            output_path="BAM/"+ ref_name + '/'), shell=True)        
                
                subprocess.call(CHIMER % dict(sample=sample, output_path="BAM/" + ref_name + '/'), shell=True)
                
                
    def variant_calling(self, bam_path, vcf_path):
        self.reference = SCRIPT_PATH + "refs/cmv_FJ527563.1_UL54.fasta"
        super().variant_calling('BAM/ul54/', 'VCF/ul54/')
        
        self.reference = SCRIPT_PATH + "refs/cmv_FJ527563.1_UL97.fasta"
        super().variant_calling('BAM/ul97/', 'VCF/ul97/')
        
    def cns(self, bam_path, cns_path, cns_x_path, min_depth_call, min_freq_thresh):
        create_dirs([cns_path+'ul54', cns_path+'ul97', cns_x_path+'ul54', cns_x_path+'ul97'])
        super().cns(bam_path+'ul54/', cns_path+'ul54/', cns_x_path+'ul54/', min_depth_call, min_freq_thresh)
        super().cns(bam_path+'ul97/',cns_path+'ul97/', cns_x_path+'ul97/', min_depth_call, min_freq_thresh)
        
    def depth(self, bam_path, depth_path):
        create_dirs([depth_path+'ul54', depth_path+'ul97/'])
        super().depth(bam_path+'ul54/', depth_path+'ul54/')
        super().depth(bam_path+'ul97/', depth_path+'ul97/')    
    
    def mafft(self, not_aligned, aligned):
        create_dirs(['alignment/ul54/', 'alignment/ul97/'])
        for ref_name, ref_path  in self.references_dict.items():
            subprocess.call(CAT % dict(path = 'CNS/'+ ref_name + '/*', 
                                       output='alignment/' + ref_name + '/all_not_aligned.fasta'), shell=True)
            utils.mafft(ref_path, 'alignment/' + ref_name + '/all_not_aligned.fasta', 
                        'alignment/' + ref_name + '/all_aligned.fasta')
        
    def qc_report(self, bam_path, depth_path, output_report):
        super().qc_report(bam_path + 'ul54/', depth_path + 'ul54/', output_report + "_" + 'ul54')
        super().qc_report(bam_path + 'ul97/', depth_path + 'ul97/', output_report + "_" + 'ul97')
        utils.create_dirs(["reports"])
        signatures.run("alignment/ul54/all_aligned.fasta", '', "reports/mutations_ul54.xlsx", show_all =  True)
        signatures.run("alignment/ul97/all_aligned.fasta", '', "reports/mutations_ul97.xlsx", show_all =  True)
        self.resistance()
    
    def resistance(self):
        ul97 = pd.read_excel('reports/mutations_ul97.xlsx')
        ul97['gene_name'] = 'UL97'
        ul54 = pd.read_excel('reports/mutations_ul54.xlsx')
        ul54['gene_name'] = 'UL54'

        mut_tbl = pd.concat([ul54,ul97]).drop(columns=['nt_position_on_genome'])

        resist = pd.read_csv(SCRIPT_PATH + "refs/CMV_GCV-R_resistance.csv").drop(columns=['nt_position_on_gene'])

        #create AA mut tbl
        filtered_columns = [col for col in mut_tbl.columns if col.endswith("_AA")]
        mut_tbl_AA = mut_tbl[['gene_name', 'aa_position_on_gene', 'R/S'] + filtered_columns].drop_duplicates()
        
        resist = pd.merge(resist, mut_tbl_AA, on=['gene_name', 'aa_position_on_gene'], how = 'left').drop_duplicates()
        
        resist.insert(0, 'has_mutation', np.where(resist['R/S'] == 'R', '+', '-'))

        #fillter mut_tbl to only contain SNPs 
        mut_tbl_snps = signatures.only_show_snp(mut_tbl)
        
        self.format_xl(mut_tbl_snps, resist)
    
    def format_xl(self, mut_tbl_snps, resist):
        # Create a Pandas Excel writer using XlsxWriter engine
        with pd.ExcelWriter('reports/mutations&resistance.xlsx', engine='xlsxwriter') as writer:
            
            
            num_samples = len(self.sample_fq_dict)
            # Write each DataFrame to a separate sheet
            resist.to_excel(writer, sheet_name='resistance_mutations', index=False)
            mut_tbl_snps.to_excel(writer, sheet_name='other_mutations', index=False)

            # Access the XlsxWriter workbook and worksheet objects
            workbook = writer.book
            worksheet1 = writer.sheets['resistance_mutations']
            worksheet2 = writer.sheets['other_mutations']

            # Define formats for conditional formatting
            gray_format = workbook.add_format({'bg_color': '#d3d3d3',
                                               'font_color': '#000000'}) 
            green_format = workbook.add_format({'bg_color': '#C6EFCE',
                                                'font_color': '#006100'})
            red_format = workbook.add_format({'font_color': 'red'})
            yellow_format = workbook.add_format({'bg_color':   '#FFFB00'})
            
            # Determine the last column index dynamically
            max_row = int(resist.shape[0])
            max_col = int(resist.shape[1])
            
            # Add conditional formatting to Sheet1
            worksheet1.conditional_format(1, 0, max_row, max_col-1, {'type': 'formula',
                                                'criteria': '=$A2="-"',
                                                'format': gray_format})
            worksheet1.conditional_format(1, 0, max_row, 0, {'type': 'formula',
                                                            'criteria': '=$A2="+"',
                                                            'format': green_format})


            max_row = int(mut_tbl_snps.shape[0])
            max_col = int(mut_tbl_snps.shape[1])
            xl_aa_start = 3+num_samples-1+2
            xl_aa_end = 3+num_samples-1+2+num_samples
            
            # Add conditional formatting to Sheet2 
            worksheet2.conditional_format(1, max_col-3, max_row, max_col-3, {'type':     'formula',
                                               'criteria': "=$" + xl_col_to_name(max_col-1) +"2=1",
                                               'format':   red_format})
            
            worksheet2.conditional_format(0,3, max_row, num_samples + 3, {'type':     'cell',
                                            'criteria': 'equal to',
                                            'value':    '"N"',
                                            'format':   gray_format})
            worksheet2.conditional_format(0,3, max_row, num_samples + 3, {'type':     'cell',
                                            'criteria': 'equal to',
                                            'value':    '"-"',
                                            'format':   gray_format})
            
            worksheet2.conditional_format(1,3, max_row, 3+num_samples-1, {'type':     'formula',
                                                'criteria': "=NOT($C2=D2)",
                                                'format':   yellow_format})
            
            worksheet2.conditional_format(0,1, max_row,  max_col, {'type':     'cell',
                                            'criteria': 'equal to',
                                            'value':    '"X"',
                                            'format':   gray_format})
            
            
            worksheet2.conditional_format(1,xl_aa_start+1, max_row, xl_aa_end, {'type':     'formula',
                                                'criteria': "=NOT($"+xl_col_to_name(xl_aa_start)+"2="+xl_col_to_name(xl_aa_start+1)+"2)",
                                                'format':   yellow_format}) 
            
            
        
        
        
        
        
        
        
        
        
        
        