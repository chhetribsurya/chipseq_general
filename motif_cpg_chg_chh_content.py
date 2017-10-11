import numpy as np
import pandas as pd
import pybedtools
import re
import os
import shutil
import time
import glob
from pyfasta import Fasta
from os.path import basename
from os.path import join
from os.path import splitext

start_time = time.time()

meme_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/MeMe/meme.txt")
fasta_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/hg_19_misc/hg19-male.fa")
fasta_idx = Fasta(fasta_file)
#query_motif = "Motif 2"

output_dir = os.path.expanduser("~/Dropbox/local_miscellaneous_data/motifs_methplot")
if not os.path.exists(output_dir):
	os.makedirs(output_dir)

#motif_output_file = "final_motifs_coordinate.bed"
final_output_file = "final_wgbs_motifs_intersect.bed"


"""Returns the list of motif_regex/motif_seq"""
def parse_meme_motif_regex(MeMe_file, motif_name=False):	
	MeMe_file = meme_file
    with open(MeMe_file, "r") as motif_file:
        meme_content = motif_file.read()
        motif_name_pattern = "Motif\s\d+\sregular\sexpression"
        motif_name_list = re.compile(motif_name_pattern).findall(meme_content)

        motif_seq_pattern = "Motif\s\d+\sregular\sexpression\n-*\n(.*?\n)"
        motif_seq_list = re.compile(motif_seq_pattern).findall(meme_content)

        motif_dict = {}
        for i in range(len(motif_name_list)):
            motif_key, motif_value = " ".join(motif_name_list[i].split()[0:2]), motif_seq_list[i].strip("\n")      
            if motif_key in motif_dict:
                motif_dict[motif_key].append(motif_seq_list)
            else:
                motif_dict[motif_key] = [motif_value]

    return(motif_dict)

motif_patterns = parse_meme_motif_regex(meme_file)
#gata_regex = motif_patterns[query_motif]


#meme_content = open(meme_file, "r").read()
def parse_meme_motif_coords(MeMe_file):
	with open(meme_file, "r") as inputfile:
		meme_content = inputfile.read()
		str_pattern  = r"Motif \d+ sites sorted by position p-value.*\n.*\n.*\n.*\n"
		motif_markers_matched_list = re.compile(str_pattern).findall(meme_content)
		header = ["chrom","start","end","strand","meme_seq", "pyfasta_seq"]
		motif_header = "\t".join(header)

		Master_dict = {}
		for each_match in motif_markers_matched_list:
			motif_key = each_match[0:8].strip()
			motif_data = []
			for i in range(len(header)):
				motif_data.append([])

			motif_match_len = len(each_match)
			start_index = meme_content.find(each_match) + motif_match_len
			end_index = meme_content.find("---", start_index)
			required_each_meme_content = meme_content[start_index: end_index]
			each_motif_lines = required_each_meme_content.strip().split("\n")
			
			for motif_line in each_motif_lines:		
				splitted = motif_line.split()
				chrom1 = splitted[0].split(":")[0]
				sequence1 = splitted[5]
				strand1 = splitted[1]
				
				if strand1 == "+":			
					start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) - 1
					end1 = start1 + len(sequence1)
					seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()
					
				if strand1 == "-":
					start1 = (int(splitted[0].split(":")[1]) + int(splitted[2])) -1 
					end1 = start1 + len(sequence1)
					seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

				req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]
				for i,item in enumerate(req):
					motif_data[i].append(item)

			motif_zip = zip(header,motif_data)
			motif_dict = dict(motif_zip)

			if motif_key not in Master_dict:
				Master_dict[motif_key] = motif_dict
		
	return(Master_dict)

Master_motif_dict = parse_meme_motif_coords(meme_file)


def combine_motif_coords_after_parsing_meme_motifs(master_motif_dict, **kwargs):
	print "motifs to be combined..."
	combine_motif_dict = {}
	count = 0
	for key, value in kwargs.iteritems():
		motif_value = "Motif "+ str(value)	
		print key, ":", motif_value
		motif_df = pd.DataFrame(master_motif_dict[motif_value])		
		count +=1			
		if motif_value not in combine_motif_dict:
			combine_motif_dict[motif_value] = motif_df
		
	concat_motif_df = pd.concat(combine_motif_dict)
	combined_motif_df = concat_motif_df.reset_index()
	final_combined_df = combined_motif_df.iloc[:,[2,6,3,7,4,5,0]]
	final_combined_df.to_csv(join(output_dir, "motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	print "\nTotal count of motifs combined:", count
	return(final_combined_df)

final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, motif_1= 1, motif_2= 2)


def alternative_combine_motif_coords_after_parsing_meme_motifs(master_motif_dict, **kwargs):
	print "motifs to be combined:\n"
	combine_motif_list = []
	count = 0
	for key, value in kwargs.iteritems():
		motif_value = "Motif "+ str(value)	
		print key, ":", original_value
		motif_df = pd.DataFrame(master_motif_dict[motif_value])
		combine_motif_list.append(motif_df.head())
		count +=1
	combined_motif_df = pd.concat(combine_motif_list, ignore_index=True)
	print "Total count of motifs combined:\n", count
	return(combined_motif_df)


### For trimming purpose:
def trim_and_parse_single_meme_motif_coords(MeMe_file, motif_number, trim_up, trim_down, trim_all_motifs=False):
	#Trim_up = 11
	#Trim_down = 0	
	if trim_all_motifs:
		Motif_num = "\d+"
	else:
		Motif_num = motif_number
	motif_value = "Motif " + str(Motif_num)
	meme_file = MeMe_file
	Trim_up = trim_up
	Trim_down = trim_down
	with open(meme_file, "r") as inputfile:
		meme_content = inputfile.read()
		str_pattern  = r"Motif " + str(Motif_num) + " sites sorted by position p-value.*\n.*\n.*\n.*\n"
		motif_markers_matched_list = re.compile(str_pattern).findall(meme_content)
		header = ["chrom","start","end","strand","meme_seq", "pyfasta_seq"]
		motif_header = "\t".join(header)

		Master_dict = {}
		for each_match in motif_markers_matched_list:
			motif_key = each_match[0:8].strip()
			motif_data = []
			for i in range(len(header)):
				motif_data.append([])

			motif_match_len = len(each_match)
			start_index = meme_content.find(each_match) + motif_match_len
			end_index = meme_content.find("---", start_index)
			required_each_meme_content = meme_content[start_index: end_index]
			each_motif_lines = required_each_meme_content.strip().split("\n")
			
			for motif_line in each_motif_lines:		
				splitted = motif_line.split()
				chrom1 = splitted[0].split(":")[0]
				strand1 = splitted[1]
				sequence1 = splitted[5]		
				Effective_len = (len(sequence1) - Trim_up)	

				if strand1 == "+":
					# fasta_idx.sequence({"chr": "chr18", "start" : (3603155 + 22 -1) + Trim_up, "stop" : (3603155 + 22 -1) + (len_3 - Trim_down), "strand" : "+"}, one_based = False).upper() # u'GATAAAG'
					start1_temp = (int(splitted[0].split(":")[1]) + int(splitted[2]) - 1)
					start1 = start1_temp + Trim_up
					end1 = start1_temp + (len(sequence1) - Trim_down)
					seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

				if strand1 == "-":
					#fasta_idx.sequence({"chr": "chr2", "start" : (112383865 + 59 -1) + Trim_down, "stop" : (112383865 + 59 -1) + (len_3 - Trim_up), "strand" : "-"}, one_based = False).upper() # u'GATAATC'
					start1_temp = ((int(splitted[0].split(":")[1]) + int(splitted[2])) -1)  
					start1 =start1_temp + Trim_down
					end1 = start1_temp + (len(sequence1) - Trim_up)
					seq_pyfasta = fasta_idx.sequence({"chr": chrom1, "start" : start1, "stop" : end1, "strand" : strand1}, one_based = False).upper()

				req = [chrom1, start1, end1, strand1, sequence1, str(seq_pyfasta)]
				for i,item in enumerate(req):
					motif_data[i].append(item)

			motif_zip = zip(header,motif_data)
			motif_dict = dict(motif_zip)
			motif_df = pd.DataFrame(motif_dict)
			if motif_key not in Master_dict:
				print motif_key
				Master_dict[motif_key] = motif_df
	#final_trimmed_df = pd.DataFrame(Master_dict[motif_value])
	final_trimmed_df = pd.concat(Master_dict).reset_index()
	return(final_trimmed_df)


df1 = trim_and_parse_single_meme_motif_coords(meme_file, 1, 11, 0)
df2 = trim_and_parse_single_meme_motif_coords(meme_file, 2, 0, 0)


def combine_trimmed_motif_df(*args):
	combine_df_list = []
	for each_df in args:
		combine_df_list.append(each_df)

	combined_trimmed_df = pd.concat(combine_df_list, ignore_index=True)
	final_trimmed_motif_df = combined_trimmed_df.iloc[:,[2,6,3,7,4,5,0]]
	final_trimmed_motif_df.to_csv(join(output_dir, "trim_based_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	return(final_trimmed_motif_df)

final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

### Using boolean style of slicing:
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("CG")]
# final_trimmed_motif_df[final_trimmed_motif_df["pyfasta_seq"].str.contains("(CG)")]

### For regex compatible:
# [m.span() for m in re.compile("GA").finditer(str_test)]
# for each in [m.span() for m in re.compile("GA").finditer(str_test)]:
# 	print each[0], each[1]
# for each_str in final_trimmed_motif_df["pyfasta_seq"].tolist():
# 	[m.span() for m in re.compile("GA").finditer(each_str)]

# df_match = final_trimmed_motif_df[~final_trimmed_motif_df["pyfasta_seq"].str.extract("(CG)").isnull()]
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([CAT]))")
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([\w+]G))").head(60)
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.contains("CG").astype(int).head(60).sum()
# df_match.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())
# df_match = final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C([A-Z]))")
# select_indices = df_match.index.get_level_values(0).tolist() ## includes repetition of index
# final_trimmed_motif_df["pyfasta_seq"].str.extractall("(C(G))")
# select_indices = CG_match_df.index.levels[0].tolist() ## includes unique indexes
# motif_seq_df = final_motif_df.iloc[select_indices]


def CG_containing_motifs_coord_df(final_motif_df_info):
	final_motif_df = final_motif_df_info
	### For finding CG containing motif dataframes, ussing boolean style of slicing:
	CG_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("CG")]
	C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]

	#CG_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([G]))").shape
	CG_regex_match_df = final_motif_df["meme_seq"].str.extractall("((CG))")
	CG_regex_grouped_df = CG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

	final_CG_motif_joined_df = pd.concat([CG_containing_motif_df, CG_regex_grouped_df], axis=1)
	final_CG_motif_joined_df = final_CG_motif_joined_df.rename(columns={0:"CG_count"}) ## since 0 based indices so add 1 to count a match:
	final_CG_motif_joined_df["CG_count"] = final_CG_motif_joined_df["CG_count"]+1 ## since 0 based indices so add 1 to count a match:
	final_CG_motif_joined_df.to_csv(join(output_dir, "CG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	
	total_CG_motif_summary_dict = {}
	total_CG_motifs = CG_containing_motif_df.shape[0]
	total_motifs = final_motif_df.shape[0]
	total_C_lacking_motifs = C_lacking_motif_df.shape[0]
	
	total_CG_motif_summary_dict["total_CG_motifs"] = total_CG_motifs
	total_CG_motif_summary_dict["total_motifs"] = total_motifs
	total_CG_motif_summary_dict["CG_motif_percent"] = (total_CG_motifs/float(total_motifs)*100)
	total_CG_motif_summary_dict["total_CG_sites_from_total_CG_motifs"] = final_CG_motif_joined_df["CG_count"].sum() #or, CG_regex_match_df.shape[0]
	total_CG_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
	total_CG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

	return(final_CG_motif_joined_df, total_CG_motif_summary_dict)

final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df)
#CG_containing_motif_df["meme_seq"].str.replace("CG","0")


"""CHG (where H is A, C or T)"""

def CHG_containing_motifs_coord_df(final_motif_df_info):

	final_motif_df = final_motif_df_info
	### For finding CHG containing motif dataframes, ussing boolean style of slicing:
	CHG_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("C[ATC]G")]
	C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]

	#CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([ATC]G))")
	CHG_regex_match_df = final_motif_df["meme_seq"].str.extractall("((C[ATC]G))")
	CHG_regex_grouped_df = CHG_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

	final_CHG_motif_joined_df = pd.concat([CHG_containing_motif_df, CHG_regex_grouped_df], axis=1)
	final_CHG_motif_joined_df = final_CHG_motif_joined_df.rename(columns={0:"CHG_count"}) ## since 0 based indices so add 1 to count a match:
	final_CHG_motif_joined_df["CHG_count"] = final_CHG_motif_joined_df["CHG_count"]+1 ## since 0 based indices so add 1 to count a match:
	final_CHG_motif_joined_df.to_csv(join(output_dir, "CHG_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	
	total_CHG_motif_summary_dict = {}
	total_CHG_motifs = CHG_containing_motif_df.shape[0]
	total_motifs = final_motif_df.shape[0]
	total_C_lacking_motifs = C_lacking_motif_df.shape[0]
	
	total_CHG_motif_summary_dict["total_CHG_motifs"] = total_CHG_motifs
	total_CHG_motif_summary_dict["total_motifs"] = total_motifs
	total_CHG_motif_summary_dict["CHG_motif_percent"] = (total_CHG_motifs/float(total_motifs)*100)
	total_CHG_motif_summary_dict["total_CHG_sites_from_total_CHG_motifs"] = final_CHG_motif_joined_df["CHG_count"].sum() #or, CHG_regex_match_df.shape[0]
	total_CHG_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
	total_CHG_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

	return(final_CHG_motif_joined_df, total_CHG_motif_summary_dict)

final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df)
#CHG_containing_motif_df["meme_seq"].str.replace("C[ATC}G","0")


"""CHH (where H is A, C or T)"""

def CHH_containing_motifs_coord_df(final_motif_df_info):

	final_motif_df = final_motif_df_info
	### For finding CHH containing motif dataframes, ussing boolean style of slicing:
	CHH_containing_motif_df = final_motif_df[final_motif_df["meme_seq"].str.contains("C[ATC][ATC]")]
	C_lacking_motif_df = final_motif_df[~final_motif_df["meme_seq"].str.contains("C")]
    C_lacking_motif_df.to_csv(join(output_dir, tf_name+"C_lacking_motif_coordinate_info.bed"), sep="\t", index=False, header=False)

	#CHH_regex_match_df = final_motif_df["meme_seq"].str.extractall("(C([ATC][ATC]))")
	CHH_regex_match_df = final_motif_df["meme_seq"].str.extractall("((C[ATC][ATC]))")
	CHH_regex_grouped_df = CHH_regex_match_df.reset_index().groupby(["level_0"]).apply( lambda x: x["match"].max())

	final_CHH_motif_joined_df = pd.concat([CHH_containing_motif_df, CHH_regex_grouped_df], axis=1)
	final_CHH_motif_joined_df = final_CHH_motif_joined_df.rename(columns={0:"CHH_count"}) ## since 0 based indices so add 1 to count a match:
	final_CHH_motif_joined_df["CHH_count"] = final_CHH_motif_joined_df["CHH_count"]+1 ## since 0 based indices so add 1 to count a match:
	final_CHH_motif_joined_df.to_csv(join(output_dir, "CHH_containing_motif_coordinate_info.bed"), sep="\t", index=False, header=False)
	total_C_lacking_motifs = C_lacking_motif_df.shape[0]
	
	total_CHH_motif_summary_dict = {}
	total_CHH_motifs = CHH_containing_motif_df.shape[0]
	total_motifs = final_motif_df.shape[0]
	total_CHH_motif_summary_dict["total_CHH_motifs"] = total_CHH_motifs
	total_CHH_motif_summary_dict["total_motifs"] = total_motifs
	total_CHH_motif_summary_dict["CHH_motif_percent"] = (total_CHH_motifs/float(total_motifs)*100)
	total_CHH_motif_summary_dict["total_CHH_sites_from_total_CHH_motifs"] = final_CHH_motif_joined_df["CHH_count"].sum() #or, CHH_regex_match_df.shape[0]
	total_CHH_motif_summary_dict["total_C_lacking_motifs"] = total_C_lacking_motifs
	total_CHH_motif_summary_dict["non_C_motif_percent"] = (total_C_lacking_motifs/float(total_motifs)*100)

	return(final_CHH_motif_joined_df, total_CHH_motif_summary_dict)

final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df)
#CHH_containing_motif_df["meme_seq"].str.replace("C[ATC][ATC]","0")



def generate_tss_binned_coords(refgene_coordinates_info, upstream_range, downstream_range, bin_size):
	#upstream = 1000
	#downstream = 1000
	#bin_size=100
	refgene_df =  refgene_coordinates_info.sort_values(["chrom","txStart","txEnd"])
	upstream = upstream_range
	downstream = downstream_range
	bin_size = bin_size
	nrows =  refgene_df.shape[0]

	bins = range(-upstream, (downstream), bin_size)
	bin_len = len(bins)
	refgene_concat_df = pd.concat([refgene_df]*bin_len, ignore_index="TRUE")
	refgene_sorted_df = refgene_concat_df.sort_values(["chrom","txStart","txEnd"])

	### Copy the bin list that is deep copy:
	bin_start_list = bins[:]
	bin_end_list = []
	for each in bin_start_list:
		bin_end_list.append(each+bin_size)

	bin_df = pd.DataFrame({"bin_start":bin_start_list, "bin_end":bin_end_list})
	bin_df = bin_df.iloc[:,[1,0]] # Arrange bin_start as first col and bin_start as 2nd col.
	bin_concat_df = pd.concat([bin_df]*nrows, ignore_index="TRUE")

	### Combine the refgene df and bin df by cbind or column wise:
	temp_refgene_df = pd.concat([refgene_sorted_df.reset_index(), bin_concat_df], axis = 1)
	temp_refgene_df["tss_midpoint"] = temp_refgene_df["txStart"]
	final_refgene_df = temp_refgene_df.loc[:,["chrom", "tss_midpoint", "bin_start","bin_end","name2","strand"]]

	""" For positive strand, i.e if strand1 == "+": 
	chrom_start = tss_midpt + (bin_start); 
	chrom_end = tss_midpt + (bin_end)"""

	refgene_pos_df = final_refgene_df.loc[final_refgene_df["strand"] == "+",]
	refgene_pos_df["chrom_start"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_start"]
	refgene_pos_df["chrom_end"] = refgene_pos_df["tss_midpoint"] + refgene_pos_df["bin_end"]
	#refgene_pos_df.iloc[:,range(1,8)]

	""" For negative strand, i.e if strand1 == "-": 
	chrom_start = tss_midpt - (bin_end ); 
	chrom_end = tss_midpt - (bin_start);
	Meth_r concept, start and end switched, so as to always maintain the higher coords for end site """

	refgene_neg_df = final_refgene_df.loc[final_refgene_df["strand"] == "-",]
	refgene_neg_df["chrom_start"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_end"]
	refgene_neg_df["chrom_end"] = refgene_neg_df["tss_midpoint"] - refgene_neg_df["bin_start"]
	#refgene_pos_df.iloc[:,range(1,8)]

	### Combine the positive and negative stranded genes:
	refgene_model_df = pd.concat([refgene_pos_df, refgene_neg_df])
	select_cols = [u'chrom', u'chrom_start', u'chrom_end', u'bin_start', u'bin_end', u'strand', u'name2', u'tss_midpoint']
	tss_coord_df1 = refgene_model_df.loc[:,select_cols]
	tss_coord_df1.to_csv(join(output_dir, "tss_coordinate_info.bed"), sep="\t", index=False, header=False)

	return(tss_coord_df1)



def load_tss_coord_pybedtool_object(file_name_with_full_path): 
	each_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["each_file"] = each_file
    sorted_peak_file = splitext(basename(each_file))[0] + "_sorted" + splitext(basename(each_file))[1]
    os.environ["sorted_peak_file"] = sorted_peak_file

    #if not os.path.exists(join(output_dir, sorted_peak_file)):
    CMD = 'sort -k1,1 -k2,2n $each_file > $output_dir/$sorted_peak_file'
    os.system(CMD)   
  
    print "Generating bedtools object..."
    peak_bed = pybedtools.BedTool(join(output_dir,sorted_peak_file))
 
    return(peak_bed)


def load_cpg_pybedtool_object(file_name_with_full_path):
	print " \nProcessing Cpg bed file\n "
	cpg_bed_file = file_name_with_full_path
	os.environ["output_dir"] = output_dir
	os.environ["cpg_bed_file"] = cpg_bed_file
	cpg_bed_edited = splitext(basename(cpg_bed_file))[0] + "_bedEdited" + splitext(cpg_bed_file)[1]
	os.environ["cpg_bed_edited"] = cpg_bed_edited

	#if not os.path.exists(join(output_dir,cpg_bed_edited)):
	CMD = '''awk 'BEGIN{FS=" "; OFS="\t"} { print $1,$2,$3+1,$5,$6,"." }' $cpg_bed_file | sort -k1,1 -k2,2n | grep -v "^*" > $output_dir/$cpg_bed_edited'''
	os.system(CMD) #subprocess.call(COMMAND, shell=True) # import subprocess
	#else:
	#    print "CpG file already converted to bed format..."

	print "Generating bedtools object..."
	Cpg_bed_file = pybedtools.BedTool(join(output_dir,cpg_bed_edited))

	return(Cpg_bed_file


def generate_tss_binned_perc_meth(refseq_final_bedfile, meth_file_list, **kwargs):
	#file_name =  kwargs["files_basename"]
	print "kwargs: ", kwargs
	refseq_bedfile = refseq_final_bedfile

	master_dict = {}

	for idx, each_file in enumerate(meth_file_list):
		file_basename = splitext(basename(each_file))[0]	
		cpg_bed_file = load_cpg_pybedtool_object(each_file)

	    """Pybedtool intersection of Cpg_bed_file and peak_bed_file"""
	    print " Processing the intersection for", file_basename
	    pybed_outfile = join(output_dir, (file_basename + "_cpg_tss_intersect.bed"))
	    pybed_outfile_v = join(output_dir, (file_basename + "_cpg_tss_not_intersect.bed"))

	    #if not os.path.exists(pybed_outfile):
	    cpg_bed_intersect = cpg_bed_file.intersect(refseq_bedfile, wa = True, wb = True, output=pybed_outfile)	    
	    print(cpg_bed_intersect.head())  
	       #Cpg_bed_file.intersect(peak_bed, wa = True, wb = True, v = True, output=pybed_outfile_v)
	    #else:
	    #	print "Pybedtools object already present"

	    ### Working with the dataframes; reading the output file of pybedtool intersect:
	    df_file = pd.read_csv(pybed_outfile, sep = "\t", header = None)
	    df_ordered = df_file.iloc[:, [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 12]]
	    df_ordered.columns = ["cpg_chrom", "cpg_start", "cpg_end", "cpg_meth", "cpg_unmeth", "peak_chrom", "peak_start", "peak_end", "bin_start", "bin_end", "gene_name" ]
	    df_ordered[["cpg_meth", "cpg_unmeth"]] = df_ordered[["cpg_meth", "cpg_unmeth"]].astype(int)
	    df_grouped =  df_ordered.groupby(["bin_start", "bin_end"]).apply(lambda x : x["cpg_meth"].sum()/float(x["cpg_meth"].sum() + x["cpg_unmeth"].sum()))
	    print "Dimension of currently intersected peak file is", df_grouped.reset_index().shape

	    if file_basename not in master_dict:
	      master_dict[file_basename] = df_grouped

		print "\nIntersection of %s completed!!!...\n\n" %(file_basename)

	### Combine all the dataframes (i.e all values of master_tf_dict) representing the tf intersects with CpG:
	cpg_intersect_combined_df = pd.concat(master_dict).reset_index()
		cpg_intersect_combined_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

	return(cpg_intersect_combined_df)



meme_file_path = "/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR*/MeMe/*.txt"
meme_file_list = glob.glob(meme_file_path)

master_tf_dict = {}
for meme_file in meme_file_list:
        #meme_file = '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL151597_SL151598_KLF6_v2[FLAG]/MeMe/meme.txt' #string instance for regex
        tf_name = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][-1] #[('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]
        SL_rep1 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][0]
        SL_rep2 = re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(meme_file)[0][1]
        print "\nCurrently processing %s tf\n" %(tf_name)
        #motif_patterns = parse_meme_motif_regex(meme_file)
        Master_motif_dict = parse_meme_motif_coords(meme_file)
        final_motif_df = combine_motif_coords_after_parsing_meme_motifs(Master_motif_dict, motif_1= 1, motif_2= 2)
        #df1 = trim_and_parse_single_meme_motif_coords(meme_file, 1, 11, 0)
        #df2 = trim_and_parse_single_meme_motif_coords(meme_file, 2, 0, 0)
        #final_trimmed_motif_df = combine_trimmed_motif_df(df1,df2)

        final_CG_containing_motif_df, final_CG_motif_summary_dict = CG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHG_containing_motif_df, final_CHG_motif_summary_dict = CHG_containing_motifs_coord_df(final_motif_df, tf_name)
        final_CHH_containing_motif_df, final_CHH_motif_summary_dict = CHH_containing_motifs_coord_df(final_motif_df, tf_name)

        if tf_name not in master_tf_dict:
                master_tf_dict[tf_name] = final_CG_motif_summary_dict
                master_tf_dict[tf_name].update(final_CHG_motif_summary_dict)
                master_tf_dict[tf_name].update(final_CHH_motif_summary_dict)

master_tf_df = pd.DataFrame(master_tf_dict)
final_master_tf_df = master_tf_df.T
final_master_tf_df.to_csv(join(output_dir, final_output_file), sep ="\t", header = True, index = True)

print "Task completed!!!"
print "\n\nTime taken to complete the analayis", time.time()-start_time



### For motif analysis using fastinterval and bedtools logic, so this is zero-based so the end element is exclusive:
bin_size = 1
upstream = 50
downstream = 50
bins = range(-upstream, (downstream+1), 1)

test_data = [] 
analysis_dict = {}


final_data = list()
final_header = ["chrom", "start", "end", "bin_start", "mid_point", "bin_end", "sequence", "motif_start", "motif_end"]

for i in range(len(final_header)):
	final_data.append(list())

coordinate_file = os.path.expanduser("~/Desktop/pub_data/python_output_files/for_gata_parse/gata_motif.bed")
with open(coordinate_file, "w") as outfile:
	print "chrom\tmid_point\tstart\tend\tbin_start\tbin_end\tsequence\tmotif_start\tmotif_end"
	for i in range(len(chrom)):
		chrome = chrom[i]
		print chrome
		motif_start = int(start[i])
		motif_end = int(end[i])
		mid_point = (motif_start + motif_end) / 2
		sequence = sequence_motif_pyfasta[i]

		for each in bins :
			bin_start = each
			bin_end = (each + bin_size)
			#For fastinterval and bedtools case:
			#bin_end = (each + bin_size)
			chrom_start = mid_point + bin_start
			chrom_end = mid_point + bin_end

			line_coords = [chrome, str(chrom_start), str(chrom_end), str(bin_start), str(mid_point), str(bin_end), sequence, str(motif_start), str(motif_end)]
			test_data.append("\t".join(line_coords))

			if chrome in analysis_dict:
				analysis_dict[chrome].append((chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords)))
			else:
				analysis_dict[chrome] = [(chrome, int(chrom_start), int(chrom_end), "\t".join(line_coords))]

				
			for i,item in enumerate(line_coords):
				final_data[i].append(item)

			print_coords = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(chrome, chrom_start, chrom_end, bin_start, mid_point, bin_end, sequence, motif_start, motif_end) 
			print print_coords
			outfile.write(print_coords + "\n")

	final_zip = zip(final_header, final_data)
	final_dict = dict(final_zip)


"""Regex match"""

In [141]: str_test
Out[141]: '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL151597_SL151598_KLF6_v2[FLAG]/MeMe/meme.txt'

In [142]: str_test1
Out[142]: '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL88833_SL88834_RAD21[FLAG]/MeMe/meme.txt'

In [143]: str_test2
Out[143]: '/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL146868_SL146869_ARID3A/MeMe/meme.txt'

In [144]: re.compile(r".*IDR_.*?_.*?_(.*?)(?=\/MeMe)").findall(str_test)
Out[144]: ['KLF6_v2[FLAG]']

In [145]: re.compile(r".*IDR_.*?_.*?_(.*?)(?=\/MeMe)").findall(str_test1)
Out[145]: ['RAD21[FLAG]']

In [146]: re.compile(r".*IDR_.*?_.*?_(.*?)(?=\/MeMe)").findall(str_test2)
Out[146]: ['ARID3A']

In [167]: re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(str_test2)
Out[167]: [('SL146868', 'SL146869', 'ARID3A')]

In [168]: re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(str_test1)
Out[168]: [('SL88833', 'SL88834', 'RAD21[FLAG]')]

In [169]: re.compile(r".*IDR_(.*?)_(.*?)_(.*?)(?=\/MeMe)").findall(str_test)
Out[169]: [('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]



In [147]: str_test4 ='/gpfs/gpfs1/home/schhetri/for_chris/batch_I/chip_analysis/IDR_SL151597_SL151598_KLF6_yuyu_v2[FLAG]/MeMe/meme/dfads/fdsfd.txt'

In [148]: re.compile(r".*IDR_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9\[\]_]+)/.*").findall(str_test4)
Out[148]: [('SL151597', 'SL151598', 'KLF6_yuyu_v2[FLAG]')]

In [149]: re.compile(r".*IDR_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9\[\]_]+)/.*").findall(str_test)
Out[149]: [('SL151597', 'SL151598', 'KLF6_v2[FLAG]')]

In [150]: re.compile(r".*IDR_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9\[\]_]+)/.*").findall(str_test1)
Out[150]: [('SL88833', 'SL88834', 'RAD21[FLAG]')]

In [151]: re.compile(r".*IDR_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9\[\]_]+)/.*").findall(str_test2)
Out[151]: [('SL146868', 'SL146869', 'ARID3A')]


In [153]: a = re.compile(r".*IDR_(?P<first>[a-zA-Z0-9]+)_(?P<second>[a-zA-Z0-9]+)_([a-zA-Z0-9\[\]_]+)\/.*").match(str_test)

In [154]: a.groups()
Out[154]: ('SL151597', 'SL151598', 'KLF6_v2[FLAG]')

In [156]: a.group('first')
Out[156]: 'SL151597'

In [157]: a.group('second')
Out[157]: 'SL151598'


In [110]: my_str = "my name is ram"

In [158]: re.compile(r"(\w+)\s+(?=ram)").findall(my_str)
Out[158]: ['is']

In [159]: re.compile(r"(?<=name)\s+(\w+)\s+(?=ram)").findall(my_str)
Out[159]: ['is']

In [160]: re.compile(r"(?<=name)\s+(\w+)").findall(my_str)
Out[160]: ['is']


# In [728]: master_dict = {}

#master_dict["tf1"] = final_CG_motif_summary_dict

# In [730]: master_dict["tf2"] = final_CHG_motif_summary_dict

# In [731]: master_dict["tf3"] = final_CHH_motif_summary_dict

# In [732]: master_dict
# Out[732]: 
# {'tf1': {'CG_motif_percent': 4.071246819338422,
#   'total_CG_motifs': 32,
#   'total_CG_sites_from_total_CG_motifs': 36,
#   'total_motifs': 786},
#  'tf2': {'CHG_motif_percent': 37.02290076335878,
#   'total_CHG_motifs': 291,
#   'total_CHG_sites_from_total_CHG_motifs': 370,
#   'total_motifs': 786},
#  'tf3': {'CHH_motif_percent': 35.75063613231552,
#   'total_CHH_motifs': 281,
#   'total_CHH_sites_from_total_CHH_motifs': 343,
#   'total_motifs': 786}}

# In [733]: pd.DataFrame(master_dict)
# Out[733]: 
#                                               tf1         tf2         tf3
# CG_motif_percent                         4.071247         NaN         NaN
# CHG_motif_percent                             NaN   37.022901         NaN
# CHH_motif_percent                             NaN         NaN   35.750636
# total_CG_motifs                         32.000000         NaN         NaN
# total_CG_sites_from_total_CG_motifs     36.000000         NaN         NaN
# total_CHG_motifs                              NaN  291.000000         NaN
# total_CHG_sites_from_total_CHG_motifs         NaN  370.000000         NaN
# total_CHH_motifs                              NaN         NaN  281.000000
# total_CHH_sites_from_total_CHH_motifs         NaN         NaN  343.000000
# total_motifs                           786.000000  786.000000  786.000000

# In [734]: pd.DataFrame(master_dict).T
# Out[734]: 
#      CG_motif_percent  CHG_motif_percent  CHH_motif_percent  total_CG_motifs  \
# tf1          4.071247                NaN                NaN             32.0   
# tf2               NaN          37.022901                NaN              NaN   
# tf3               NaN                NaN          35.750636              NaN   

#      total_CG_sites_from_total_CG_motifs  total_CHG_motifs  \
# tf1                                 36.0               NaN   
# tf2                                  NaN             291.0   
# tf3                                  NaN               NaN   

#      total_CHG_sites_from_total_CHG_motifs  total_CHH_motifs  \
# tf1                                    NaN               NaN   
# tf2                                  370.0               NaN   
# tf3                                    NaN             281.0   

#      total_CHH_sites_from_total_CHH_motifs  total_motifs  
# tf1                                    NaN         786.0  
# tf2                                    NaN         786.0  
# tf3                                  343.0         786.0  






### Dictionary updates:
# In [771]: master_tf_dict["tf_3"] =  final_CG_motif_summary_dict

# In [772]: master_tf_dict["tf_3"].update(final_CHH_motif_summary_dict)

# In [773]: master_tf_dict["tf_3"].update(final_CHG_motif_summary_dict)

# In [774]: master_tf_dict
# Out[774]: 
# {'tf_1': {'CG_motif_percent': 4.071246819338422,
#   'CHG_motif_percent': 37.02290076335878,
#   'CHH_motif_percent': 35.75063613231552,
#   'total_CG_motifs': 32,
#   'total_CG_sites_from_total_CG_motifs': 36,
#   'total_CHG_motifs': 291,
#   'total_CHG_sites_from_total_CHG_motifs': 370,
#   'total_CHH_motifs': 281,
#   'total_CHH_sites_from_total_CHH_motifs': 343,
#   'total_motifs': 786},
#  'tf_2': {'CG_motif_percent': 4.071246819338422,
#   'CHG_motif_percent': 37.02290076335878,
#   'CHH_motif_percent': 35.75063613231552,
#   'total_CG_motifs': 32,
#   'total_CG_sites_from_total_CG_motifs': 36,
#   'total_CHG_motifs': 291,
#   'total_CHG_sites_from_total_CHG_motifs': 370,
#   'total_CHH_motifs': 281,
#   'total_CHH_sites_from_total_CHH_motifs': 343,
#   'total_motifs': 786},
#  'tf_3': {'CG_motif_percent': 4.071246819338422,
#   'CHG_motif_percent': 37.02290076335878,
#   'CHH_motif_percent': 35.75063613231552,
#   'total_CG_motifs': 32,
#   'total_CG_sites_from_total_CG_motifs': 36,
#   'total_CHG_motifs': 291,
#   'total_CHG_sites_from_total_CHG_motifs': 370,
#   'total_CHH_motifs': 281,
#   'total_CHH_sites_from_total_CHH_motifs': 343,
#   'total_motifs': 786}}

# In [775]: pd.DataFrame(master_tf_dict)
# Out[775]: 
#                                              tf_1        tf_2        tf_3
# CG_motif_percent                         4.071247    4.071247    4.071247
# CHG_motif_percent                       37.022901   37.022901   37.022901
# CHH_motif_percent                       35.750636   35.750636   35.750636
# total_CG_motifs                         32.000000   32.000000   32.000000
# total_CG_sites_from_total_CG_motifs     36.000000   36.000000   36.000000
# total_CHG_motifs                       291.000000  291.000000  291.000000
# total_CHG_sites_from_total_CHG_motifs  370.000000  370.000000  370.000000
# total_CHH_motifs                       281.000000  281.000000  281.000000
# total_CHH_sites_from_total_CHH_motifs  343.000000  343.000000  343.000000
# total_motifs                           786.000000  786.000000  786.000000

# In [776]: pd.DataFrame(master_tf_dict).T
# Out[776]: 
#       CG_motif_percent  CHG_motif_percent  CHH_motif_percent  total_CG_motifs  \
# tf_1          4.071247          37.022901          35.750636             32.0   
# tf_2          4.071247          37.022901          35.750636             32.0   
# tf_3          4.071247          37.022901          35.750636             32.0   

#       total_CG_sites_from_total_CG_motifs  total_CHG_motifs  \
# tf_1                                 36.0             291.0   
# tf_2                                 36.0             291.0   
# tf_3                                 36.0             291.0   

#       total_CHG_sites_from_total_CHG_motifs  total_CHH_motifs  \
# tf_1                                  370.0             281.0   
# tf_2                                  370.0             281.0   
# tf_3                                  370.0             281.0   

#       total_CHH_sites_from_total_CHH_motifs  total_motifs  
# tf_1                                  343.0         786.0  
# tf_2                                  343.0         786.0  
# tf_3                                  343.0         786.0  


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result
This function will work in Python 2 and 3 for all dicts. e.g. given dicts a to g:

z = merge_dicts(a, b, c, d, e, f, g) 
and key value pairs in g will take precedence over dicts a to f, and so on.



