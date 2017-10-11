import pandas as pd
import os
import re
from scipy import stats
from pybedtools import BedTool
from glob import glob
from os.path import join
from os.path import splitext
from os.path import basename

input_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/split_chrom_intersect_hotspot_gerp_scores")
output_dir = join(input_dir, "analysis_output")

if not os.path.exists(output_dir):
	os.makedirs(output_dir)


chromwise_gerp_analysis_files = glob(join(input_dir, "*intersect.bed"))
chromwise_gerp_analysis_files_df = []

for each in chromwise_gerp_analysis_files:
	df = pd.read_csv(each, sep="\t", header=None)
	df_groupby = df.groupby([0,1,2,3,4])[10].apply(lambda x: x.mean())
	df_groupby = df_groupby.reset_index()
	df_groupby.columns = ["chrom", "start", "end", "anno_state", "tf_count", "mean_rs_score"]
	
	bins = [1,2,3,4,5,10,20,30,40,50,70,100,500]
	names = ["1", "2", "3", "4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
	df_groupby['binned_tf_count'] = pd.cut(df_groupby["tf_count"], bins, right=False, labels=names)
	df_groupby_sorted = df_groupby.sort_values(["tf_count"])
	chromwise_gerp_analysis_files_df.append(df_groupby_sorted)

combined_gerp_df = pd.concat(chromwise_gerp_analysis_files_df, ignore_index=True)
combined_gerp_df.to_csv(join(output_dir, "combined_chromwise_gerp_analysis_boxplot_data.txt"), sep="\t", header=True, index=False)


### Gerp elements analysis:
output_dir = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/split_chrom_intersect_hotspot_gerp_elements")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

gerp_element_file = os.path.expanduser("~/for_chris/batch_I/gerp_analysis/split_chrom_intersect_hotspot_gerp_elements/combined_hotspot_gerp_elements_intersect.bed")
gerp_elem_df = pd.read_csv(gerp_element_file, sep="\t", header=None)
bins = [1,5,10,20,30,40,50,70,100,500]
names = ["1-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-69", "70-99", "100+"]
gerp_elem_df['binned_tf_count'] = pd.cut(gerp_elem_df.iloc[:,4], bins, right=False, labels=names)
gerp_elem_df_sorted = gerp_elem_df.sort_values([4]).reset_index(drop=True)
gerp_elem_df_sorted.to_csv(join(output_dir, "combined_gerp_element_analysis_boxplot_data.txt"), sep="\t", header=True, index=False)


""" Perform 2 sample KS test to see if the distribution b/w 2 samples is significantly different:"""
input_file = os.path.expanduser("~/Dropbox/encode_3/gerp_analysis/analysis_output/combined_gerp_element_analysis_boxplot_data_with_gerp_elem_overlap_fraction_final_edited.txt")
gerp_df = pd.read_csv(input_file, sep="\t")

""" Boundary of seperation for low and highly bound regions """
tf_number = 80

gerp_df_lowly_bound = gerp_df.loc[gerp_df["TF_bound"] <= tf_number]
gerp_df_highly_bound = gerp_df.loc[~(gerp_df["TF_bound"] <= tf_number)]

gerp_df.shape[0] == gerp_df_lowly_bound.shape[0] + gerp_df_highly_bound.shape[0]
#gerp_lowly_bound_array = gerp_df_lowly_bound["TF_bound"].as_matrix()
#gerp_highly_bound_array = gerp_df_highly_bound["TF_bound"].as_matrix()


""" KS (Kolmogorov-Smirnov) 2 sample test for promoter associated """
gerp_df_lowly_bound_promoter = gerp_df_lowly_bound.loc[gerp_df_lowly_bound["re_anno"] == "Promoter_assoc"]
gerp_df_highly_bound_promoter = gerp_df_highly_bound.loc[gerp_df_highly_bound["re_anno"] == "Promoter_assoc"]

gerp_df_lowly_bound_promoter_array = gerp_df_lowly_bound_promoter["fraction_overlap_score_per_100bp"].as_matrix()
gerp_df_highly_bound_promoter_array = gerp_df_highly_bound_promoter["fraction_overlap_score_per_100bp"].as_matrix()

### KS test:
ks_2_sample_test  = stats.ks_2samp(gerp_df_lowly_bound_promoter_array, gerp_df_highly_bound_promoter_array)
ks_pval_prom = ks_2_sample_test[1]

if ks_pval_prom == 0:
	new_ks_pval_prom = 1e-323
	print "KS pvalue for promoter assoc :", new_ks_pval_prom
else:
	new_ks_pval_prom = ks_pval_prom
	print "KS pvalue for promoter assoc :", new_ks_pval_prom



""" KS (Kolmogorov-Smirnov) 2 sample test for enhancer associated """
gerp_df_lowly_bound_enhancer = gerp_df_lowly_bound.loc[gerp_df_lowly_bound["re_anno"] == "Enhancer_assoc"]
gerp_df_highly_bound_enhancer = gerp_df_highly_bound.loc[gerp_df_highly_bound["re_anno"] == "Enhancer_assoc"]

gerp_df_lowly_bound_enhancer_array = gerp_df_lowly_bound_enhancer["fraction_overlap_score_per_100bp"].as_matrix()
gerp_df_highly_bound_enhancer_array = gerp_df_highly_bound_enhancer["fraction_overlap_score_per_100bp"].as_matrix()

### KS test:
ks_2_sample_test  = stats.ks_2samp(gerp_df_lowly_bound_enhancer_array, gerp_df_highly_bound_enhancer_array)
ks_pval_enh = ks_2_sample_test[1]

if ks_pval_enh == 0:
	new_ks_pval_enh = 1e-323
	print "KS pvalue for enhancer assoc :", new_ks_pval_enh
else:
	new_ks_pval_enh = ks_pval_enh
	print "KS pvalue for enhancer assoc :", new_ks_pval_enh

### Final check if all the data points were included:
gerp_df_lowly_bound.shape[0] == gerp_df_lowly_bound_promoter.shape[0] + gerp_df_lowly_bound_enhancer.shape[0]
gerp_df_highly_bound.shape[0] == gerp_df_highly_bound_promoter.shape[0] + gerp_df_highly_bound_enhancer.shape[0]
gerp_df.shape[0] == gerp_df_highly_bound_promoter.shape[0] + gerp_df_lowly_bound_promoter.shape[0] + gerp_df_highly_bound_enhancer.shape[0] + gerp_df_lowly_bound_enhancer.shape[0] 




