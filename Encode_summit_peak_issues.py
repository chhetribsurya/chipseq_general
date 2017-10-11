import glob
import pandas as pd
import os
from os.path import join
from os.path import splitext
from os.path import basename

file_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
file_list = glob.glob(file_dir)

concat_df_list = []
for each_file in file_list:
    df = pd.read_csv(each_file, sep="\t", header=None)
    # df_cols = df.iloc[0:2,] # for just 2 rows
    # df_cols = df.iloc[0:(df.shape[0]+1),]
    # df_cols = df
    df_final = df.iloc[:, [0,1,2,9]]
    df_final.columns = ["chrom", "start", "end", "summit_peak"]
    
    # mimic the summit peak from the peak file:
    df_final["peak_center"] = (df_final["end"] - df_final["start"])/2
    concat_df_list.append(df_final)
    print df_final

combined_df = pd.concat(concat_df_list)
combined_df.shape[0] == len(combined_df["summit_peak"].tolist()) # check if the length are equal
#combined_df.shape[0] == combined_df["summit_peak"].shape[0] # check if the length are equal

# check if the center of peak is exactly the same as peak summit:
combined_df["summit_peak"] == combined_df["peak_center"]
if combined_df["summit_peak"].tolist() == combined_df["peak_center"].tolist():
	print "\nGreat! The peak center has the exact coordinate as the peak summit"
else:
	print "\nSorry! The peak center and the peak summit doesn't match"

