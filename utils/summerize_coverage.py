#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 07:53:18 2022

@author: hagar
"""

def summerize_coverage(report_path):
    output = report_path.split(".csv")[0]+"_coverage_read_count.xlsx"
    
    import pandas as pd
    
    report = pd.read_csv(report_path).sort_values(by=['sample'],ignore_index=True)
    
    coverage_df = pd.DataFrame(columns=["sample", "Sabin1", "Sabin2", "Sabin3", "nOPV2"])
    read_count_df = pd.DataFrame(columns=["sample","Sabin1", "Sabin2", "Sabin3", "nOPV2"])
    coverage_col_name = "coverage%" if "contig_based" in report_path else "coverage_CNS_5%"
    i=0
    while i < (len(report)-1):
        sample = report.iloc[i]["sample"].split("-")[0].split(".")[0]
        s1_cv = ""
        s1_rc = ""
        s2_cv = ""
        s2_rc = ""
        s3_cv = ""
        s3_rc = ""
        opv_cv= ""
        opv_rc = ""
        if "19.1" in report.iloc[i]["sample"] or "Sabin1" in report.iloc[i]["sample"]  or "S1" in report.iloc[i]["sample"]:
            s1_cv = report.iloc[i][coverage_col_name]
            s1_rc = report.iloc[i]["mapped_reads"]
            i = i+1 if i < (len(report)-1) else len(report)-1
        if "20.1" in report.iloc[i]["sample"] or "Sabin2" in report.iloc[i]["sample"] or "S2" in report.iloc[i]["sample"]:
            s2_cv = report.iloc[i][coverage_col_name]
            s2_rc = report.iloc[i]["mapped_reads"]
            i = i+1 if i < (len(report)-1) else len(report)-1
        if "21.1" in report.iloc[i]["sample"] or "Sabin3" in report.iloc[i]["sample"] or "S3" in report.iloc[i]["sample"]:
            s3_cv = report.iloc[i][coverage_col_name]
            s3_rc = report.iloc[i]["mapped_reads"]
            i = i+1 if i < (len(report)-1) else len(report)-1
        if "nOPV2" in report.iloc[i]["sample"]:
            opv_cv = report.iloc[i][coverage_col_name]
            opv_rc = report.iloc[i]["mapped_reads"]
            i = i+1 if i < (len(report)-1) else len(report)-1
        coverage_df=coverage_df.append({"sample":sample, "nOPV2": opv_cv,"Sabin1":s1_cv, "Sabin2":s2_cv, "Sabin3":s3_cv},ignore_index=True)
        read_count_df=read_count_df.append({"sample":sample, "nOPV2": opv_rc,"Sabin1":s1_rc, "Sabin2":s2_rc, "Sabin3":s3_rc},ignore_index=True)
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
        coverage_df.to_excel(writer, sheet_name='coverage', index = False)
        read_count_df.to_excel(writer, sheet_name='read_count', index = False)
        writer.save()
