#!/usr/bin/env python

import pandas as pd
import numpy as np

#df = pd.read_excel("~/Dropbox/HepG2_prom_enh_ctcf_muliple_associated_tfs.xls", sep="\t")
df = pd.read_excel("~/Dropbox/HepG2_prom_enh_ctcf_muliple_associated_tfs_final.xlsx", sep="\t")
#df_filtered = df.dropna(subset=["Ideas Prom", "Ideas Enh", "Ideas CTCF", "Ideas Other", "IDEAS multiple"], how="all")
df_long = pd.melt(df.loc[:,list(df.columns[0:1]) + list(df.columns[5:])], id_vars=["all_TFs_HepG2"])
df_final = df_long[pd.notnull(df_long["value"])]
df_final.columns = ["TF_name", "annotation", "value"]
df_final = df_final[~df_final["TF_name"].duplicated()]


# excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xls"
# df_xls = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")
# df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
# xls_tf_list =  df_xls["Target"].tolist()

df_xls = pd.read_excel("~/Dropbox/xls_files/Encode_full_hepg2_datasets.xlsx", sep="\t", sheetname="unique_TFs_DBF_CR_annotation")
df_xls["target_tf_edit"] = df_xls["Target"].apply(lambda x : str(x).split("_")[0].split("[FLAG]")[0])
df_xls_final = df_xls.loc[:,["Target", "Category", "target_tf_edit"]]


### Merge the 2 dataframes:
merged_df = pd.merge(df_xls_final, df_final, left_on="target_tf_edit", right_on="TF_name", how="left")
df_select = merged_df.loc[:,["Target", "Category", "annotation"]]
df_select["annotation"] = df_select["annotation"].fillna("Unknown")
#df_select["annotation"].replace({np.nan:"Unknown"})
df_select["annotation"].replace({"IDEAS multiple": "Ideas multiple"}, inplace=True)
df_select["label"] =  df_select["annotation"].replace({'Ideas Enh':"orange", 'Ideas Prom':"red", 'Ideas Other':"green", 'Ideas multiple': "burlywood", 'Unknown':"azure", 'Ideas CTCF':"blueviolet"})
df_select.to_csv("~/Dropbox/TFs_Annotation_file.txt", sep="\t", header=True, index=False)
