#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 11:02:05 2023

@author: hagar
"""

import pandas as pd

path = "/mnt/project1/projects/HIV/seq_27072023/QC/"

fasta = pd.read_csv(path + "fasta.csv")
forma = pd.read_excel(path + "format.xlsx")

mergi = pd.merge(forma, fasta, on=["SAMPLE_No_NGS","Region_Protein"], how = "left")

mergi.to_excel(path + "merged.xlsx", index=False)
