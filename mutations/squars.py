#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 14:22:20 2022

@author: hagar
"""

import sys
import pandas as pd
from utils.utils import get_sequences, hamming_distance
import seaborn as sns
import matplotlib.pyplot  as plt

#given sequences return a df of differences- 1=>different 0=identical


def dist_mat(sequences):
    distmat= pd.DataFrame(columns=sequences.keys(), index=sequences.keys())
    for sample1, record1 in sequences.items():
        for sample2, record2 in sequences.items():
            distmat[sample1][sample2]=hamming_distance(sequences[sample1],sequences[sample2])
    return distmat.astype('int32')

def squars(distmat, title, names):
    plt.clf()
    plt.figure(figsize=(25, 20))
    res=sns.heatmap(distmat, cmap="YlGnBu", yticklabels=names, xticklabels=names)
    res.set_yticklabels(res.get_ymajorticklabels(), fontsize = 6)
    res.set_xticklabels(res.get_xmajorticklabels(), fontsize = 6)
    plt.title(title, fontsize=16)
    plt.savefig("images/"+title+".png")

def mutations(sequences, start, end, names, title):
        differences= pd.DataFrame(columns=sequences.keys(), index=range(end-start))
        for sample1, record1 in sequences.items():
            for i in range(end-start):
                differences[sample1][i] = 0
                for sample2, record2 in sequences.items():
                    if not record1[i] == record2[i]:
                        differences[sample1][i] += 1 
        differences = differences.replace(1, 0)
                        
        y_axis = []
        for name in names:
            #temp comment to remove precents
            # sum = differences[name].count()
            # precents=(sum/(end-start)*100)
            # y_axis.append(name+": "+str(int(precents))+ "%")
            y_axis.append(name)
        differences = differences.astype('int32')
        plt.figure(figsize=(8, 3))
        res = sns.heatmap(differences.T, cmap="YlGnBu", cbar=False, yticklabels=y_axis)
        res.set_xticklabels(res.get_xmajorticklabels(), fontsize=5.5)
        res.set_yticklabels(res.get_ymajorticklabels(), fontsize=8)
        plt.title(title, fontsize=16)
        plt.savefig("images/"+title + " mutations" + ".png")


if __name__ == "__main__":
    alignment_path = sys.argv[1]
    title = sys.argv[2]
    ref = False
    
    sequences = get_sequences(alignment_path)
    
    for sample, seq in sequences.items():
        sequences[sample] = list(seq)
    
    dist_mut = dist_mat(sequences)
    squars(dist_mut, title + " Hamming Distance", sequences.keys())
    seq_length = len(next(iter(sequences.values())))
    mutations(sequences, 0, seq_length, sequences.keys(),title + " mutations")


#to sort the distance mut:
# for col in data.columns:
#     data = data.sort_values(by = col)
    
# final = pd.DataFrame()

# for index in data.index:
#     final[index] = data[index]

