#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:03:49 2022

@author: hagar
"""

from utils.utils import get_sequences
import pandas as pd
'''
not in use
'''

fasta = "/mnt/project1/projects/POLIO/VDPV3_2022/fasta/sabin1/all_aligned.fasta"
seqs = get_sequences(fasta)

#fix polio name
seqs_new = {}
for sample, seq in seqs.items():
    if not "Sabin1" in sample:
        print(sample)
        new_name = sample.split("_")[1]
        if len(new_name) > 2:
            new_name = new_name.lower().split("polio")[1]
        if len(new_name) > 2:
            new_name = new_name.split(".")[0]
        seqs_new[new_name] = list(seq)
        print(new_name)

df = pd.DataFrame(seqs_new)
ref = pd.Series(list(seqs["Sabin1"]))
df.columns = df.columns.str.replace("_NT", "")
df.columns = df.columns.astype(int)

df = df.sort_index(axis=1)

mutation_count = []
last_seq = ref
for column in df:
    current_seq = df[column]
    compare = current_seq.compare(last_seq)
    #drop "N" and "-"
    compare = compare.drop(compare[compare.self == "-"].index)
    compare = compare.drop(compare[compare.other == "-"].index)
    
    compare = compare.drop(compare[compare.self == "N"].index)
    compare = compare.drop(compare[compare.other == "N"].index)
    
    mutation_count.append((compare))
    last_seq = current_seq


