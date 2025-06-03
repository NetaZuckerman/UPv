#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 05:53:52 2023

@author: hagar
"""
import argparse


def parser():
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', help='fastq folder')
        parser.add_argument('-r','--reference' ,help='reference')
        parser.add_argument('-t','--threads',default='12' ,help='set max threads for a prosses. default = 12')
        parser.add_argument('-mr','--multi_ref', action='store_true', help='multi fasta reference') #
        parser.add_argument('--drop_joint_reads', action='store_true', help='drop reads mapped to more than one refenerece when mapping against multi-fasta reference. default: False') #
        parser.add_argument('-f', '--flu', action='store_true', help="influenza segements analysis") #store_true will store the argument as true
        parser.add_argument('--de_novo', nargs='?', const=True, default=False, help='Perform de-novo analysis')
        parser.add_argument('--polio', action='store_true', help="PolioVirus analysis") #store_true will store the argument as true
        parser.add_argument('--HIV', help="HIV analysis for genes: PR, RT and Integrase. do not provide reference sequnce. provide format table.") #store_true will store the argument as true
        parser.add_argument('--sensitive', action='store_true', help="sensitive local alignment for HIV pipeline.") #store_true will store the argument as true
        parser.add_argument('--HSV', action='store_true', help="Herpes 1 analysis for genes: UL23, UL30, UL42. do not provide reference sequnce.") #store_true will store the argument as true
        parser.add_argument('--CMV', action='store_true', help="Cetomegalovirus (human herpesvirus 5) analysis") #store_true will store the argument as true
        parser.add_argument('--covid', action='store_true', help="covid19 variants analysis") #store_true will store the argument as true
        parser.add_argument('--HAV', action='store_true', help="Hepatitis A analysis. do not provide reference sequnce.") #store_true will store the argument as true
        parser.add_argument('--minion', help="use minimap2 to analyse minion reads. provide barcodes.csv (csv with 2 columns: barcode|sample)") #store_true will store the argument as true
        parser.add_argument('-gb','--gb_file' ,help='insert gb file and get reference regions report')
        parser.add_argument('-m','--mutations_table' , action='store_true',help='mutations table reprort. gb file flag is mandatory')
        parser.add_argument('--annotations' ,help='add annotations.csv to merge with the mutation table')
        parser.add_argument('-rg','--regions_file' ,help='insert gene regions file for mini (csv with 4 columns: GENE | START | END | STRAND)')
        parser.add_argument('--mini' , action='store_true',help='run only mutation analysis. this flag requires --input flag as alignment file.')
        parser.add_argument('-v','--vcf' , action='store_true',help='generates vcf files.')
        parser.add_argument('--cns_min_freq_thresh' ,default='0.6', help='Minimum frequency threshold(0 - 1) to call consensus. (Default: 0.6)') 
        parser.add_argument('--cns_min_depth_call' ,default='5', help='Minimum depth to call consensus. (Default: 5). Note- CNS folder will contain minimum depth of 1. CNS_X will contain minimum depth of X.') 
        args = parser.parse_args()
        return args.reference, args.input, args.threads, args.flu, args.de_novo, args.polio, args.CMV, args.HIV, \
                   args.sensitive, args.HSV, args.minion, args.gb_file, args.regions_file, args.mutations_table, \
                      args.annotations, args.mini, args.vcf,args.cns_min_freq_thresh, args.cns_min_depth_call, \
                           args.multi_ref, args.drop_joint_reads
