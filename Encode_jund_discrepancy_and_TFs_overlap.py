#!/usr/bin/env python

import pandas as pd
import glob
import os
import subprocess
import shutil
from distutils.dir_util import copy_tree
import re
from os.path import join
from os.path import basename
from collections import defaultdict


output_dir = os.path.expanduser("~/for_chris/TF_comparison_ab_flag_and_diff_labs")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

excel_file = os.path.expanduser("~/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx")
# df_xls = pd.read_excel(excel_file, sep="\t", sheet_name="Sheet2")

""" All TF containing dir """
#dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
all_tf_file_list = glob.glob(dir_path) # represents all tf file list

dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/dups/SL*"
dup_file_list = glob.glob(dir_path) # represents all tf file list

"""Read multiple excel sheets from the same file"""
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx"

dup_xls_df = pd.ExcelFile(excelfile).parse("dups_DBF_CR_annotation")
#dup_xls_df["Category_split_join"] = dup_xls_df["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
dup_tf_list =  dup_xls_df["Target"].dropna().tolist()

dup_dir =

hepg2_file_list = [ each_file for each_file in glob.glob(join(ideas_file_dir,"*_ideas_whole_genome")) 
				if os.path.basename(each_file).startswith("HepG2")]

all_tf_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total"
tf_file_list = glob.glob(join(all_tf_dir,"SL*"))


all_tf_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total"


master_tf_dict = {}
for each in dup_tf_list:
	splitted_tf =  re.split("[ _ []", each)[0]
	if splitted_tf not in master_tf_dict:
		master_tf_dict[splitted_tf] = glob.glob(join(all_tf_dir, "SL*" + splitted_tf + "*"))


del master_tf_dict["RAD21"][2]

del master_tf_dict["CTCF"][1]

del master_tf_dict["POLR2A"][1:5]

del master_tf_dict["CEBPB"][1]















uniq_xls_df = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")
#df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
unique_tf_list =  uniq_xls_df["Target"].dropna().tolist()


""" Check the xls TF list and the file suffix tf_list to make sure, if they are in the TF list"""
TF_check_list = []
for each_file in all_tf_file_list:
	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
	TF_check_list.append(tf_name)

# df_xls[df_xls["Target"].isin(TF_check_list)] # or,
if sorted(TF_check_list) == sorted(xls_tf_list):
	print "\nGreat!! TF list resembles b/w the file list and xls tf list..."
else:
	print "\nWarning: Your files TF list doesn't resemble the xls tf list..."


""" Select the category of TF for chosing the file list """
cr_tf_file_list = []
for each_tf in cr_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			cr_tf_file_list.append(each_file)

dbf_tf_file_list = []
for each_tf in dbf_df["Target"].tolist():
	for each_file in all_tf_file_list:
		if each_file.endswith(each_tf):
			dbf_tf_file_list.append(each_file)



output_dir = os.path.expanduser("~/for_chris/batch_I/idr_passed_peaks_total")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

sub_dir_dup = join(output_dir, "dups")
if not os.path.exists(sub_dir_dup):
	os.makedirs(sub_dir_dup)


sub_dir_flag = join(output_dir, "flagged")
if not os.path.exists(sub_dir_flag):
	os.makedirs(sub_dir_flag)


sub_dir_unique = join(output_dir, "unique_TFs")
if not os.path.exists(sub_dir_unique):
	os.makedirs(sub_dir_unique)


meme_output_dir = os.path.expanduser("~/for_chris/batch_I/chip_motif_analysis_total")
if not os.path.exists(meme_output_dir):
    os.makedirs(meme_output_dir)

sub_meme_dir_flagged = join(meme_output_dir, "flagged")
if not os.path.exists(sub_meme_dir_flagged):
	os.makedirs(sub_meme_dir_flagged)

sub_meme_dir_dup = join(meme_output_dir, "dups")
if not os.path.exists(sub_meme_dir_dup):
	os.makedirs(sub_meme_dir_dup)

sub_meme_dir_unique = join(meme_output_dir, "unique_TFs")
if not os.path.exists(sub_meme_dir_unique):
	os.makedirs(sub_meme_dir_unique)

### The excel file which is processed with duplicates and flagged details:
### Output of script peak_count_of_peak_files_from_excel_TFlist.py
excel_file = os.path.expanduser("~/for_chris/Encode_full_hepg2_datasets.xlsx")

"""Read the duplicate file from multiple excel sheets of the same file"""
df_xls = pd.ExcelFile(excel_file).parse("dups")
df_xls["Target"] = df_xls["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls = df_xls.dropna()

### cutting down the df to make it manageable:
df_xls_select = df_xls.iloc[:,0:6]


""" Duplicate files and dir transfer"""

### Transfer the duplicate TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0]; #print file_name
	shutil.move(file_name, sub_dir_dup)


### Transfer the duplicate motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	shutil.move(dir_name, sub_meme_dir_dup)




######################
######################
######################



"""Read the flagged file from multiple excel sheets of the same file"""
df_xls_flagged = pd.ExcelFile(excel_file).parse("flagged")
df_xls_flagged["Target"] = df_xls_flagged["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls_flagged = df_xls_flagged.dropna()

### cutting down the df to make it manageable:
df_xls_flagged_select = df_xls_flagged.iloc[:,0:6]

""" Flagged files and dir transfer"""

### Transfer the flagged TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_flagged_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0] ; #print file_name
	shutil.move(file_name, sub_dir_flag)


### Transfer the flagged motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_flagged_select.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	shutil.move(dir_name, sub_meme_dir_flagged)



#########################
#########################
#########################

def copyDirectory(src, dest):
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
        print('Directory not copied. Error: %s' % e)


"""Read the unique_TF file from multiple excel sheets of the same file"""
df_xls = pd.ExcelFile(excel_file).parse("unique_TFs")
df_xls["Target"] = df_xls["Target"].apply(lambda x: "_".join(str(x).split("-"))) 
df_xls = df_xls.dropna() # df_xls = df_xls.dropna(how="all")

### cutting down the df to make it manageable:
df_xls_unique = df_xls.iloc[:,0:6]


""" unique files and dir transfer"""

### Transfer the unique TF files using row by rows TF_name and SL#:
peak_file_path = "~/for_chris/batch_I/idr_passed_peaks_total"

for index,row in df_xls_unique.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(peak_file_path, SL_rep1 + "*" + SL_rep2 + "*" ))
	file_name = glob.glob(full_path)[0]; #print file_name
	shutil.copy(file_name, sub_dir_unique)


### Transfer the unique motif dir:
motif_file_path = "~/for_chris/batch_I/chip_motif_analysis_total"

for index,row in df_xls_unique.iterrows():
	SL_rep1 = row["Exp ID rep1"]
	SL_rep2 = row["Exp ID rep2"]
	tf_name = row["Target"]
	lab_name = row["Lab"]
	print SL_rep1, SL_rep2, tf_name, lab_name
	full_path = os.path.expanduser(join(motif_file_path, "IDR" + "*" + SL_rep1 + "*" + SL_rep2 + "*" ))
	dir_name = glob.glob(full_path)[0]; #print dir_name
	
	os.environ["dir_name"] = dir_name
	os.environ["output_dir"] = sub_meme_dir_unique
	os.system("cp -r $dir_name $output_dir")

	



	
