#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 06:12:04 2022

@author: hagar
"""


#%% 
import pandas as pd
from datetime import datetime

months = {"1" : "JAN",
         "2" : "FEB",
         "3" : "MAR",
         "4" : "APR",
         "5" : "MAY",
         "6" : "JUN",
         "7" : "JUL",
         "8" : "AUG",
         "9" : "SEP",
         "01" : "JAN",
         "02" : "FEB",
         "03" : "MAR",
         "04" : "APR",
         "05" : "MAY",
         "06" : "JUN",
         "07" : "JUL",
         "08" : "AUG",
         "09" : "SEP",
         "10" : "OCT",
         "11" : "NOV",
         "12" : "DEC"}
log = pd.read_csv("/home/hagar/polio_log.csv")

column_names = ["strain", "date", "location"]
df = pd.DataFrame(columns = column_names)

headers = []

for index, row in log.iterrows():
    row["Collection date"] = str(row["Collection date"])
    Date = row["Collection date"]
    if row["Collection date"] and not row["Collection date"] == "?" and not row["Collection date"] == "nan" and "." in row["Collection date"]:
        date = row["Collection date"].split(".")
        
        if(len(date[2]) == 2):
            Date = datetime.strptime(row["Collection date"], '%d.%m.%y').date()
        else:
            Date = datetime.strptime(row["Collection date"], '%d.%m.%Y').date()
        
        row["Collection date"] = str(date[0]) + months[str(date[1])] + str(date[2])
    header = "ISR_" + row["ngs_number"] + "_" + str(row["Sample ID"]) + "_" + str(row["Collection date"]) + "_" + str(row["Place"]).replace("-","_")
    strain = header.replace(" ", "_")
    
    location = str(row['Place']).replace('Jerusalem','JER')
    df.loc[len(df.index)] = [strain, Date, location]
    
    headers.append(header.replace(" ", "_"))
    
    
    
#%%
from utils import get_sequences

output_path = "/home/hagar/"
london = get_sequences("/mnt/project1/projects/POLIO/VDPV2_2022/run10/wg_analysis/map_to_london/alignment/all_not_aligned.fasta")
s2 = get_sequences("/mnt/project1/projects/POLIO/VDPV2_2022/run10/wg_analysis/map_to_s2/alignment/all_not_aligned.fasta")
P1 = get_sequences("/mnt/project1/projects/POLIO/VDPV2_2022/run10/P1_analysis/alignment/all_not_aligned.fasta")
seqs_22 = {}
seqs_21 = {}
p1_seqs_22 = {}
p1_seqs_21 = {}
#change headers 

for sample, seq in london.items():
    sample_name = sample.replace("Consensus_","")
    for header in headers:
        if sample_name in header:
            date = header.split("_")[3]
            if date.endswith("22"):
                seqs_22[header] = london[sample]
                p1_seqs_22[header] = P1[sample]
            elif date.endswith("21"):
                seqs_21[header] = s2[sample]
                p1_seqs_21[header] = P1[sample]
            else:
                print(sample + " date issue")
            break


def write_fasta(seqs, output_file): 
    with open(output_file, 'w') as f:
        for header, seq in seqs.items():
            f.write(">" + header + '\n')
            f.write(seq + '\n')
        
write_fasta(seqs_21, output_path + "run10_2021.fasta")
write_fasta(seqs_22, output_path + "run10_2022.fasta")
write_fasta(p1_seqs_21, output_path + "run10_P1_2021.fasta")
write_fasta(p1_seqs_22, output_path + "run10_P1_2022.fasta")


