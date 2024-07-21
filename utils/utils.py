#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 11:03:13 2022

@author: hagar
"""
import os
import sh
import subprocess
import shutil
import multiprocessing as mp
from Bio import SeqIO 
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import csv 

MAFFT = os.path.dirname(__file__)+"/MAFFT.sh %(not_aligned)s %(reference)s %(aligned)s"
RM_SPADES = "find spades -type f ! -name 'transcripts.fasta' -delete" # removes spades output expet of the contigs file 
#split bam file
SPLIT = "bamtools split -in %(bam)s -reference"
ambiguous_nucleotides = ["W", "Y", "R", "S", "D","K","M","V","H","B","X"]


def get_sample_fq_dict(fastq_path):
    '''
    generate dict of sample short name and its fastq (R1) path.

    Parameters
    ----------
    fastq_path : str
        path to fastq folder.

    Returns
    -------
    sample_fq : dict {sample : R1 fastq path}

    '''
    sample_fq = {}
    skip_files=["R2", "Undetermined", "unpaired", "singletons"]
    for r1 in os.listdir(fastq_path):
        if any(skip_file in r1 for skip_file in skip_files) or "fast" not in r1:
            continue
        sample = r1.split("_")[0].split(".fastq")[0] #sample short name
        #raise error if the sample name contains 'R1' 
        if 'R1' in sample:
            raise ValueError("Sample name should not contain 'R1'.\nCheck your fastq files names.")
        sample_fq[sample] = r1
    return sample_fq 

 
def split_bam(dir):
    '''
    split bam files in dir by reference.

    Parameters
    ----------
    dir : str
        path to bam folder.

    Returns
    -------
    None.

    '''
    for bam_file in os.listdir(dir):        
        if "sorted" in bam_file and "bai" not in bam_file:
            subprocess.call(SPLIT % dict(bam=dir+bam_file), shell=True)
            os.remove(dir + bam_file)
  
    
def create_dirs(dirs):
    '''
    create directories. if folder exist, remove and recreate.

    Parameters
    ----------
    dirs : list
        list of directories to create.

    Returns
    -------
    None.

    '''
    for dir in dirs:
        if os.path.exists(dir):
            shutil.rmtree(dir,ignore_errors=True)
        os.makedirs(dir,exist_ok=True)
  
          
def remove_from_name(dir, to_remove):
    '''
    remove str from files name in a givan directory.

    Parameters
    ----------
    dir : str
        the directory of the files.
    to_remove : str
        the string to remove from the file names.

    Returns
    -------
    None.

    '''
    for file in os.listdir(dir):
        file = dir + file
        os.rename(file, file.replace(to_remove, ''))


def change_header(dir):
    '''
    change first header of fasta files in a given directory to the file name.

    Parameters
    ----------
    dir : str
        directory of fasta files.

    Returns
    -------
    None.

    '''
    for file in os.listdir(dir):
        new_header = ">" + file.split(".fa")[0]
        file = dir + file        
        sh.sed("-i", "1s/.*/" + new_header + "/", file)
        

def mafft(reference,not_aligned, aligned):
    '''
    preform multiple alignemnt using augur align. by running a bash script "MAFFT.sh".

    Parameters
    ----------
    reference : str
        path to reference sequence (all sequences will be aligned to the reference).
    not_aligned : str
        path to not aligned multi-fasta file.
    aligned : str
        path to aligned multi-fasta file - the output.

    Returns
    -------
    None.

    '''
    subprocess.call(MAFFT % dict(not_aligned=not_aligned, reference=reference, aligned=aligned), shell=True)


def rm_spades():
    '''
    remove all spades output exept of the contigs fasta.

    '''
    subprocess.call(RM_SPADES, shell=True)


def get_sequences(alignment_file):
    '''
    read multi-fasta file and save it in a dictionary

    Parameters
    ----------
    alignment_file : str
        path to fasta file.

    Returns
    -------
    sequences : dict {header : sequence}

    '''
    sequences = {}
    alignment = SeqIO.to_dict(SeqIO.parse(alignment_file, 'fasta'))
    for sample, record in alignment.items():
        sequences[sample] = str(record.seq).upper()

        sequences[sample.replace("Consensus_", "").split("_threshold")[0]] = sequences.pop(sample)

    return sequences


def mutations_positions(sequences, no_n = 1):
    '''
    get mutations positions list by comparing all sequences to each other.

    Parameters
    ----------
    sequences : dict {sample : sequence}
    no_n : BOOL, optional
        ignore N's when =1. The default is 1.

    Returns
    -------
    mutations_positons : list
        list of position where mutation.

    '''
    mutations_positons = []
    seq_length = len(next(iter(sequences.values())))
    for pos in range(seq_length-1):
        temp = ""
        for sample, record in sequences.items():
            if not temp:
                temp = record[pos]
            if no_n and record[pos] in ["N"]:
                break
            if record[pos] in ambiguous_nucleotides:
                break
            if not temp == record[pos] :
                mutations_positons.append(pos+1)
                break
    return mutations_positons

def run_mp(threads, func, arg):
    with mp.Pool(threads) as pool:
        pool.map(func,arg)
        pool.close()
        pool.join()
 
def hamming_distance(seq1, seq2):
    '''
    calculate hamming distance of 2 sequences.
    ignoring N's and gaps

    Parameters
    ----------
    seq1 : str
        sequence.
    seq2 : str
        sequence.

    '''
    df = pd.DataFrame()
    df["seq1"] = pd.Series(seq1)
    df["seq2"] = pd.Series(seq2)
    df['difference'] = np.where((df["seq1"] == df["seq2"]) | (df["seq1"] == "N") | (df["seq2"] == "N") | (df["seq1"] == "-") | (df["seq2"] == "-"), 0, 1)
    
    return (df["difference"].sum())


def write_sub_fasta(fasta, path, regions, gene, strand='+'):
    '''
    write part of a sequence as fasta in a given path

    Parameters
    ----------
    fasta : str
        genomic sequence.
    path : str
        path to outpus file.
    regions : tuple
        gene regions (start, end).
    gene : str
        gene name.
    strand : char, optional
        '-' for complement stand. The default is '+'.

    Raises
    ------
    ValueError
        strand must be + or -.

    '''
    with open(path + gene + ".fasta",'w') as f:
        start = regions[0]
        end = regions[1]
        for header, seq in fasta.items():
            if strand=='+':
                sub_seq = seq[start-1:end-1]
            elif strand=='-':
                sub_seq = str(Seq(seq[start:end]).reverse_complement())
            else:
                raise ValueError("strand must be + or -.")
            f.write(">" + header + '\n')
            f.write(sub_seq + '\n')    
    

def fix_cns_header(path):
    '''
    iterates all fasta files in a given path and simplify their headers (removes spaces).

    Parameters
    ----------
    path : str
        fasta path.


    '''
    for file in os.listdir(path):
        if file.endswith(".fasta") or file.endswith(".fa"):
            fasta = get_sequences(path + file)
            
            ofile = open(path + file, "w")
            
            for header, seq in fasta.items():
                ofile.write(">" + header + "\n" + seq + "\n")
    
            ofile.close()
            
def get_barcodes(barcode_csv):
    '''
    get barcode|sample table for minion run from csv.

    Parameters
    ----------
    barcode_csv : str
        path to barcodes.csv file.

    barcodes : dict
        {sample: barcode}.

    '''
    with open(barcode_csv, mode='r') as infile:
        reader = csv.reader(infile)
        barcodes = dict((rows[1],rows[0]) for rows in reader)
    barcodes.pop('sample')
    return barcodes

translate_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}
