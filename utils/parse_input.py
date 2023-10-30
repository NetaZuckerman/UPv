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
        parser.add_argument('-t','--treads',default='12' ,help='set max treads for a prosses. default = 12')
        parser.add_argument('-mr','--multi_ref', action='store_true', help='multi fasta reference') #
        parser.add_argument('-f', '--flu', action='store_true', help="influenza segements analysis") #store_true will store the argument as true
        parser.add_argument('-d', '--de_novo', action='store_true', help="de-novo analysis") #store_true will store the argument as true
        parser.add_argument('--polio', action='store_true', help="PolioVirus analysis") #store_true will store the argument as true
        parser.add_argument('--HIV', help="HIV analysis for genes: PR, RT and Integrase. do not provide reference sequnce. provide format table.") #store_true will store the argument as true
        parser.add_argument('--HSV', action='store_true', help="HSV analysis for genes: UL23, UL30, UL42. do not provide reference sequnce.") #store_true will store the argument as true
        parser.add_argument('-c', '--cmv', action='store_true', help="cetomegalovirus (human herpesvirus 5) analysis") #store_true will store the argument as true
        parser.add_argument('-gb','--gb_file' ,help='insert gb file and get reference regions report')
        parser.add_argument('-m','--mutations_table' , action='store_true',help='mutations table reprort. gb file flag is mandatory')
        parser.add_argument('-rg','--regions_file' ,help='insert gene regions file for mini')
        parser.add_argument('--mini' , action='store_true',help='run only mutation analysis. this flag requires --input flag as alignment file.')
        parser.add_argument('--skip_spades' , action='store_true',help='skip spades analysis used in polio and de novo classes. turn on this flag only if you already run spades once')
        parser.add_argument('-v','--vcf' , action='store_true',help='generates vcf files.')
        parser.add_argument('--cnsThresh' ,default='0.6', help='Minimum frequency threshold(0 - 1) to call consensus. (Default: 0.6)') #defualt is set in the main script (upy.py)
        args = parser.parse_args()
        return args.reference, args.input, args.treads, args.flu, args.de_novo, args.polio, args.cmv, args.HIV, \
                   args.HSV, args.gb_file, args.regions_file, args.mutations_table, \
                       args.mini, args.skip_spades, args.vcf, args.cnsThresh, args.multi_ref
        