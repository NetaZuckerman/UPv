#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 14:22:20 2022

@author: hagar
"""

import sys
import pandas as pd
from utils.utils import get_sequences, hamming_distance
import numpy as np
import seaborn as sns
import matplotlib.pyplot  as plt
#given sequences return a df of differences- 1=>different 0=identical

def seq_difference(seq1, seq2):
 df = pd.DataFrame()
 df["seq1"] = pd.Series(seq1)
 df["seq2"] = pd.Series(seq2)
 df['difference'] = np.where((df["seq1"] == df["seq2"]) | (df["seq1"] == "N") | (df["seq2"] == "N") | (df["seq1"] == "-") | (df["seq2"] == "-"), 0, 1)
 
 return df['difference']

def dist_mat(sequences):
    distmat= pd.DataFrame(columns=sequences.keys(), index=sequences.keys())
    for sample1, record1 in sequences.items():
        for sample2, record2 in sequences.items():
            distmat[sample1][sample2]=hamming_distance(sequences[sample1],sequences[sample2])
    return distmat.astype('int32')

def squars(distmat, title, names):
    plt. clf()
    res=sns.heatmap(distmat, cmap="YlGnBu", yticklabels=names, xticklabels=names)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 8.5)
    
    plt.title(title, fontsize=16)
    plt.savefig("images/"+title+".png")

def mutations(sequences, start, end, names, title):
        differences= pd.DataFrame(columns=["CMV_1","CMV_2"], index=range(end-start))
        differences["CMV_1"] = seq_difference(sequences["CMV_1"], sequences["herpesvirus5"])
        differences["CMV_2"] = seq_difference(sequences["CMV_2"], sequences["herpesvirus5"])
        
        # for i in range(end-start):
        #     differences['CMV_1'][i] = 0
        #     differences['CMV_2'][i] = 0
            
            # if not sequences['CMV_1'][i] == sequences['herpesvirus5'][i]:
            #     differences['CMV_1'][i] = 1
            # if not sequences['CMV_2'][i] == sequences['herpesvirus5'][i]:
            #     differences['CMV_2'][i] = 1
        

        
                        
        y_axis = []
        for name in names:
            sum = differences[name].sum()
            precents=(sum/(end-start)*100)
            y_axis.append(name+": "+str(round(precents,2))+ "%")
        differences = differences.astype('int32')
        plt.figure(figsize=(30, 4))
        res = sns.heatmap(differences.T, cmap="YlGnBu", cbar=False, yticklabels=y_axis)
        res.set_xticklabels(res.get_xmajorticklabels(), fontsize=7)
        res.set_yticklabels(res.get_ymajorticklabels(), fontsize=10)
        plt.title(title, fontsize=16)
        plt.savefig("images/"+title + " differences" + ".png")
        #differences.to_csv("images/"+title + " differences" + ".csv")


if __name__ == "__main__":
    alignment_path = sys.argv[1]
    title = sys.argv[2]
    sequences = get_sequences(alignment_path)
    
    
    
    for sample, seq in sequences.items():
        sequences[sample] = list(seq)
    
    # dist_mut = dist_mat(sequences)
    # squars(dist_mut, "Hamming Distance", sequences.keys())
    seq_length = len(next(iter(sequences.values())))
    mutations(sequences, 0, seq_length, ["CMV_1","CMV_2"],"mutations")
    





