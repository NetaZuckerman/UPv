#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:50:48 2023

@author: hagar
"""

ev_path = "/data3/code/references/Enteroviruses/Enteroviruses.fasta"
qc_path = "/mnt/project1/projects/POLIO/PanEV/pilot2_09112023/analysis/QC/QC_report.csv"

import pandas as pd

pd.read_csv(qc_path)

with open(ev_path , "r") as f:
    evs = f.readlines()
evs = [x for x in evs if x.startswith(">")]

df = pd.DataFrame(columns=["reference_id", "reference_name", "enterovirus_type"])
for ev in evs:
    ref_id = ev.split(' ')[0][1:]
    ev_type = ev.strip()[-1]
    ref_name = ''
    names = ev.split(' ')
    for i in range(len(names)):
        if names[i].lower().startswith(('cox', 'echo', 'display', 'species')):
            if names[i].lower().startswith(('display', 'species')):
                ref_name = names[i].split("=")[1]
            else:
                ref_name = names[i] + names[i+1]
    if not ref_name:
        for i in range(len(names)):
            if names[i].lower().startswith("entero"):
                ref_name = names[i].strip() + names[i+1].strip()
    
    df.loc[len(df.index)] = [ref_id, ref_name, ev_type] 

df.to_csv(ev_path.replace("fasta","csv"), index=False)