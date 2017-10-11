import pandas as pd
import numpy as np
import pybedtools
import pickle
import scipy 
import os
import re
import glob
import scipy.stats as stats
import statsmodels.stats.multitest as smm
from rpy2.robjects.packages import importr  #stats = importr('stats')
from rpy2.robjects.vectors import FloatVector
from os.path import join
from os.path import basename
from os.path import splitext


null_sequence_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/unique_hepg2_analysis_total"
main_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/luciferase_enhancer_assay_analysis"
#main_dir = os.path.expanduser("~/Dropbox/for_chris/batch_I/hepg2_liver_gene_sp_analysis")
output_dir =  join(main_dir, "files")

if not os.path.exists(main_dir):
	os.makedirs(main_dir)
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

#refgene_file = "/gpfs/gpfs1/home/schhetri/for_chris/refGene_hg19"
#refgene_file = os.path.expanduser("~/Dropbox/for_chris/refGene_hg19")

""" All TF containing dir """
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
#dir_path = os.path.expanduser("~/Dropbox/local_miscellaneous_data/chip_hepg2_tf/idr_passed_peaks_total/unique_TFs/SL*")
all_tf_file_list = glob.glob(dir_path) # represents all tf file list


""" Read hepg2 luciferase enhancer elements genes """
excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/HepG2_luciferase_enhancer_assay.xlsx"
#excelfile = os.path.expanduser("~/Dropbox/for_chris/Liver_specific_genes_HepG2.xlsx")
df_xls = pd.ExcelFile(excelfile).parse()

""" Extract the coordinates for luciferase test elements """
luciferase_coord_df = df_xls.loc[:,["chr", "start", "end", "activity"]]
luciferase_activity_list = df_xls["activity"].unique().tolist()

#for each in luciferase_activity_list:
#	print luciferase_coord_df[luciferase_coord_df["activity"] == each]


# """Read multiple excel sheets from the same file"""
# excelfile = "/gpfs/gpfs1/home/schhetri/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx"
# #excelfile = os.path.expanduser("~/Dropbox/for_chris/Encode_full_hepg2_datasets_DBF_CR.xlsx")
# df_xls = pd.ExcelFile(excelfile).parse("unique_TFs_DBF_CR_annotation")
# df_xls["Category_split_join"] = df_xls["Category"].apply(lambda x: "_".join(str(x).split("/"))) 
# xls_tf_list =  df_xls["Target"].tolist()

# cr_df = df_xls[df_xls["Category"] == "CR/CF"]
# dbf_df = df_xls[~(df_xls["Category"] == "CR/CF")]

# """ Check the xls TF list and the file suffix tf_list to make sure, if they are in the TF list"""
# TF_check_list = []
# for each_file in all_tf_file_list:
# 	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
# 	TF_check_list.append(tf_name)

# if sorted(TF_check_list) == sorted(xls_tf_list):
# 	print "\nGreat!! TF list resembles b/w the file list and xls tf list..."
# else:
# 	print "\nWarning: Your files TF list doesn't resemble the xls tf list..."

# """ Select the category of TF for chosing the file list """
# cr_tf_file_list = []
# for each_tf in cr_df["Target"].tolist():
# 	for each_file in all_tf_file_list:
# 		if each_file.endswith(each_tf):
# 			cr_tf_file_list.append(each_file)


# dbf_tf_file_list = []
# for each_tf in dbf_df["Target"].tolist():
# 	for each_file in all_tf_file_list:
# 		if each_file.endswith(each_tf):
# 			dbf_tf_file_list.append(each_file)


# #######################
# """ Parse the ideas state"""

# # ideas_hepg2_file = os.path.expanduser("/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state/HepG2_ideas_whole_genome")
# ideas_mnemonics_file = os.path.expanduser("~/for_encode/spp/analysis_scripts/ideas_table.txt")

# """ Select the state number, to be analysed on """
# #select_state_num = range(1,9+1) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:
# select_state_num = range(1,9) + [15,16] + range(17,20+1) # Add +1; range() is exclusive:

# """ Generate file list """
# ideas_file_dir = "/gpfs/gpfs1/home/schhetri/encode_datasets/ideas_state"


# output_dir = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/unique_hepg2_analysis_total"
# if not os.path.exists(output_dir):
# 	os.makedirs(output_dir)


# """ Read the mnemonics file, and create a dictionary of states with state number ID"""
# read_file = pd.read_csv(ideas_mnemonics_file, sep="\t")
# state_list = read_file["Mnemonics"]
# state_num_list = [ i for i in range(1,len(state_list)+1) ]
# ideas_state_dict = dict(zip(state_num_list,state_list))
# target_state = [ideas_state_dict[each] for each in select_state_num] #Depends on the select_state_num variable input

#######################


""" Preparing for significant peak and background peaks overlap """
# dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks/test_analysis/SL*"
dir_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs/SL*"
tf_file_list = glob.glob(dir_path)

null_generate_script = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_generate.py"
null_parameters = "-x 1 -r 1 "
null_hg19_indices = "/gpfs/gpfs1/home/schhetri/kmers_svm/surya_script/nullseq_indices_hg19/"

master_luciferase_dict = {}
master_sigcount_dict = {}
master_bgcount_dict = {}
total_peak_count_dict  = {}
total_luciferase_dict = {}


for each_file in tf_file_list:
	tf_luciferase_dict = {}
	tf_name = re.compile(r".*narrowPeak_(.*)$").findall(basename(each_file))[0]
	print "\nCurrently Processing %s\n" %(tf_name)
	df_read = pd.read_csv(each_file, sep="\t", header=None)
	df_read = df_read.iloc[:,[0,1,2]]
	df_read.columns = ["tf_chrom", "start", "end"]
	df_read["mid_point"] = (df_read["start"] + df_read["end"])/2
	df_read["tf_start"] = df_read["mid_point"].astype(int) - (50)
	df_read["tf_end"] = df_read["mid_point"].astype(int) + (50)
	df_read["tf_name"] = tf_name
	total_peak_count_dict[tf_name] = df_read.shape[0]


	""" Preparing the file for significant peak overlap"""
	select_cols = ["tf_chrom","tf_start","tf_end","tf_name"]
	df_read = df_read.loc[:,select_cols]
	tf_bedfile = pybedtools.BedTool.from_dataframe(df_read)
	print "Preparing the intersection for", tf_name 
	print df_read


	""" Preparing the file for background peak overlap; bg stands for background"""
	#select_cols_bg = ["tf_chrom","tf_start","tf_end"]
	#df_bg_read = df_read.loc[:,select_cols_bg]
	#print df_bg_read
	df_read.to_csv(join(null_sequence_dir, "pre_background_file_" + tf_name), sep="\t", header=False, index=False)
	os.environ["null_generate_script"] = null_generate_script
	os.environ["null_parameters"] = null_parameters
	os.environ["null_hg19_indices"] = null_hg19_indices
	os.environ["input_file"] = join(null_sequence_dir, "pre_background_file_"+tf_name)
	os.environ["output_file"] = join(null_sequence_dir, "background_file_"+tf_name)
	

	if not os.path.exists(join(null_sequence_dir, "background_file_"+tf_name)):
		print "Generating the null sequence for", tf_name
		os.system("python $null_generate_script $null_parameters -o $output_file $input_file hg19 $null_hg19_indices")
	else:
		print "Background peak file exists for", tf_name
	#print "Preparing the background intersection for", tf_name 
	tf_bg_bedfile = pybedtools.BedTool(join(null_sequence_dir, "background_file_"+tf_name))
	print tf_bg_bedfile.head()

	
	""" Intersection b/w significant and each states, and likewise for backgound peaks """ 
	sig_count_list = []
	sig_state_list = []
	bg_count_list = []
	bg_state_list = []

	for each_state in luciferase_activity_list:
		each_uniq_df = luciferase_coord_df[luciferase_coord_df["activity"].isin([each_state])]
		each_uniq_bedfile = pybedtools.BedTool.from_dataframe(each_uniq_df)
		print each_uniq_bedfile.head()

		if each_state not in total_luciferase_dict:
			total_luciferase_dict[each_state] = each_uniq_df.shape[0]

		""" For significant peak overlap"""
		print "\nCurrently processing the intersection for.... %s and %s" %(tf_name, each_state)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_state + "_unique_ideas_intersect.bed"))
		tf_gene_intersect = tf_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_gene_intersect.count() > 0:
			each_intersect_df = pd.read_csv(tf_gene_intersect.fn, sep="\t", header=None)
			print tf_gene_intersect.head()
		else:
			each_intersect_df = None

		### For the first time, need to create a dictionary, the update/merge the new dictionary into this
		### first time created dict: tf_luciferase_dict[each_state] = {"sig_intersect": each_intersect_df}
		tf_luciferase_dict[each_state] = {"sig_intersect": each_intersect_df}
		tf_luciferase_dict[each_state].update({"sig_count": tf_gene_intersect.count()})
		sig_count_list.append(tf_gene_intersect.count())
		sig_state_list.append(each_state)


		""" For background peak overlap; bg stands for background"""		
		print "\nCurrently processing the background intersection for.... %s and %s" %(tf_name, each_state)
		pybed_outfile = join(output_dir, (tf_name + "_" + each_state + "_background_intersect.bed"))
		tf_bg_gene_intersect = tf_bg_bedfile.intersect(each_uniq_bedfile, wa=True, wb=True, output=pybed_outfile)

		if tf_bg_gene_intersect.count() > 0:
			each_bg_intersect_df = pd.read_csv(tf_bg_gene_intersect.fn, sep="\t", header=None)
			print tf_bg_gene_intersect.head()
		else:
			each_bg_intersect_df = None
					
		tf_luciferase_dict[each_state].update({"bg_intersect": each_bg_intersect_df})
		tf_luciferase_dict[each_state].update({"bg_count": tf_bg_gene_intersect.count()})
		bg_count_list.append(tf_bg_gene_intersect.count())
		bg_state_list.append(each_state)

	master_sigcount_dict[tf_name] = dict(zip(sig_state_list, sig_count_list))
	master_bgcount_dict[tf_name] = dict(zip(bg_state_list, bg_count_list))
	master_luciferase_dict[tf_name] = tf_luciferase_dict

with open(join(output_dir,"final_master_luciferase_dict.pkl"), "w") as outfile:
	pickle.dump(master_luciferase_dict, outfile)

### To load back the pickled python objects:
# with open(join(output_dir,"final_master_luciferase_dictas_dict_3tf.pkl")) as infile:
#	df = pickle.load(infile)
master_sig_bg_dict = {}
for each_tf in master_sigcount_dict:
	sig_key_list =  master_sigcount_dict[each_tf].keys()
	sig_value_list = master_sigcount_dict[each_tf].values()
	bg_key_list =  master_bgcount_dict[each_tf].keys()
	bg_value_list = master_bgcount_dict[each_tf].values()
	total_peak_count = total_peak_count_dict[each_tf]
	total_activity_list = total_luciferase_dict.keys()
	total_sites_list = total_luciferase_dict.values()
	master_sig_bg_dict[each_tf]=pd.DataFrame({"sig_state":sig_key_list, "sig_hits":sig_value_list, "bg_state": bg_key_list, "bg_hits": bg_value_list, \
												 "total_peaks":total_peak_count, "activity":total_activity_list, "total_activity_sites": total_luciferase_dict.values()})

combine_tf_df = pd.concat(master_sig_bg_dict)	
final_tf_df = combine_tf_df.reset_index().drop("level_1", axis=1)	
final_tf_df = final_tf_df.rename(columns={"level_0" : "tf_name"})
select_cols = ["tf_name", "activity", "sig_hits", "bg_hits", "total_activity_sites", "total_peaks"]
final_ideas_df = final_tf_df.loc[:,select_cols]
final_ideas_df.to_csv(join(output_dir,"final_tf_ideas_peak_count_with_bg.bed"), sep="\t", header=True, index=False)

with open(join(output_dir,"final_luciferase_df.pkl"), "w") as outfile:
	pickle.dump(final_ideas_df, outfile)

with open(join(output_dir,"final_luciferase_df.pkl")) as infile:
    final_ideas_df = pickle.load(infile)

final_ideas_df["sig_hits_unbound"] = final_ideas_df["total_peaks"] - final_ideas_df["sig_hits"]
final_ideas_df["bg_hits_unbound"] = final_ideas_df["total_peaks"] - final_ideas_df["bg_hits"]

# calculate fisher exact test & bh/bonferroni correction for multiple hypothesis/test correction:
def fishers_test(df_rows):
	x1 = df_rows["sig_hits"]
	x2 = df_rows["sig_hits_unbound"]
	y1 = df_rows["bg_hits"]
	y2 = df_rows["bg_hits_unbound"]
	pvalue = stats.fisher_exact([[x1,y1],[x2,y2]], "greater")[1] # R: 
	# pvalue = fisher.test(rbind(c(1,9),c(11,3)), alternative="greater")$p.value
	# by default the R's fisher exact has "lesser" as alt. hypothesis. So, to match
	# with python, set the alt hypothesis as; alternative = "less"

	""" #total_pval_test = df_rows["sig_hits"].shape[0]
		#bonferroni_correction = min(pvalue*total_pval_test, 1)
		#return(pvalue, df_rows["tf_name"]) """

	return pd.Series({"pvalue": pvalue, "sig_hits": x1, "bg_hits": y1, "tf_name": df_rows["tf_name"], 
						"activity": df_rows["activity"], "total_activity_sites": df_rows["total_activity_sites"], "total_peaks": df_rows["total_peaks"] })

final_ideas_pval_df = final_ideas_df.apply(fishers_test, axis=1)

""" Some scipy.stats statistical insights with pearson correlation, and Kolmogorov Smirnov"""
""" # df_test= final_ideas_df.apply(fishers_test, axis=1)
# final_df = df_test.apply(pd.Series), or; final_df = pd.DataFrame(df_test.values.tolist())
# df_clean = df.dropna()
# stats.pearsonr(df_clean['VIXCLS'], df_clean['GDP'])
# (0.7205766921228921, 0.48775429164459994)
# Where the first value in the tuple is the correlation value, and second is the p-value
# 2 sample Kolmogorov Smirnov test,
# scipy.stats.ks_2samp(dataset1, dataset2) 
# scipy.stats.chisquare([11105, 13573], f_exp=[698402, 918898]) """

# Replace the enrichment value 0 with the highest decimal place the python can generate i.e 1.0e-323
final_ideas_pval_df = final_ideas_pval_df.sort_values(["pvalue"])
final_ideas_pval_df["pvalue"] = final_ideas_pval_df["pvalue"].replace("0", 1.0e-323)

pvalue_list = final_ideas_pval_df["pvalue"].tolist()
rej, pval_corr = smm.multipletests(pvalue_list, alpha=0.05, method="fdr_bh")[:2] # method="bonferroni" or "hs"; if needed
final_ideas_pval_df["pval_corr"] = list(pval_corr)

order_cols = ["tf_name", "activity", "sig_hits", "bg_hits", "pvalue", "pval_corr", "total_activity_sites", "total_peaks"]
final_pval_qval_df =  final_ideas_pval_df.loc[:,order_cols]
#final_pval_qval_df["percent_overlap"] = (final_pval_qval_df["sig_hits"]/final_pval_qval_df["total_peaks"])*100
#sig_hits_sorted_df = final_pval_qval_df.sort_values("sig_hits", ascending=False)

# Log transformation:
final_pval_qval_df["-log10(pvalue)"] = -np.log10(final_pval_qval_df["pvalue"].astype(float))
final_pval_qval_df["-log10(pval_corr)"] = -np.log10(final_pval_qval_df["pval_corr"].astype(float))
final_pval_qval_df_sorted = final_pval_qval_df.replace(-0, 0)

with open(join(output_dir,"final_pval_qval_df.pkl"), "w") as outfile:
	pickle.dump(final_pval_qval_df_sorted, outfile)

with open(join(output_dir,"final_pval_qval_df.pkl")) as infile:
    final_pval_qval_df_sorted = pickle.load(infile)

# """ Enrichment value of 0.001 """
# final_pval_qval_001 = final_pval_qval_df_sorted[final_pval_qval_df_sorted["pval_corr"] < 0.001]
# final_pval_qval_001.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.001.bed"), sep="\t", header=True, index=False)
# final_pval_qval_heatmap_001 = final_pval_qval_001.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
# final_pval_qval_heatmap_001.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.001.bed"), sep="\t", header=True, index=True)

# # Based on the output of heatmap decided to exclude some of the datas:
# final_pval_qval_state_filtered_001 = final_pval_qval_001[final_pval_qval_001["states"].isin(["Enh", "EnhF", "PromCtcf"])]
# final_pval_qval_state_filtered_001.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_0.001_state_filtered.bed"), sep="\t", header=True, index=True)
# final_pval_qval_heatmap_001 = final_pval_qval_state_filtered_001.pivot(index="tf_name", columns="states", values = "-log10(pval_corr)")
# final_pval_qval_heatmap_001.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_0.001_state_filtered.bed"), sep="\t", header=True, index=True)

""" Enrichmnent value heatmap """
final_pval_qval_df_sorted.to_csv(join(output_dir,"final_tf_luciferase_enh_enrichment_with_pval_qval.bed"), sep="\t", header=True, index=False)
final_pval_qval_heatmap_01 = final_pval_qval_df_sorted.pivot(index="tf_name", columns="activity", values = "-log10(pval_corr)")
final_pval_qval_heatmap_01.to_csv(join(output_dir,"final_tf_luciferase_enh_enrichment_pval_corr_heatmap_data.bed"), sep="\t", header=True, index=True)

# # Based on the output of heatmap decided to exclude some of the datas:
# sig_gene_list = final_pval_qval_df[final_pval_qval_df["pvalue"] < 0.05]["activity"].unique().tolist()
# final_pval_qval_gene_filtered = final_pval_qval_df_sorted[final_pval_qval_df_sorted["activity"].isin(sig_gene_list)]
# final_pval_qval_gene_filtered.to_csv(join(output_dir,"final_tf_enrichment_with_pval_qval_gene_filtered.bed"), sep="\t", header=True, index=True)
# final_pval_qval_heatmap = final_pval_qval_gene_filtered.pivot(index="tf_name", columns="activity", values = "-log10(pval_corr)")
# final_pval_qval_heatmap.to_csv(join(output_dir,"final_tf_enrichment_pval_corr_heatmap_data_gene_filtered.bed"), sep="\t", header=True, index=True)





