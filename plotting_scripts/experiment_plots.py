import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
from scipy.stats import kde
import sys

from RT import *
from plots import *

def find_targets_decoys_on_threshold(data_df, feature_name, round_flag, keep_nth_element_n):
	if round_flag==1:
		data_df[feature_name] = data_df[feature_name].round(4)
	
	unique_scores = data_df[feature_name].unique()
	unique_scores = np.sort(unique_scores)
	
	unique_scores = unique_scores[::keep_nth_element_n]
	
	num_targets = []
	num_decoysPerTargets = []
	print(len(unique_scores))
	counter = 0
	for score in unique_scores:
		if counter % 100 == 0:
			print(counter, end=', ')
		targets_N = float(len(data_df.loc[(data_df[feature_name] <= score) & (data_df['Decoy'] == 1)][feature_name].tolist()))
		decoys_N = float(len(data_df.loc[(data_df[feature_name] <= score) & (data_df['Decoy'] == -1)][feature_name].tolist()))
		num_targets.append(targets_N)
		if targets_N == 0:
			num_decoysPerTargets.append(0.0)
		else:
			num_decoysPerTargets.append(decoys_N / targets_N)
			if decoys_N / targets_N >= 0.05:
				break
		counter += 1
		
	return [num_targets, num_decoysPerTargets]	
	
	
def evalue_number_of_hits_databases(PSMs_extended_df_1, PSMs_ensembl_ref_df_1, PSMs_uniprot_iso_df_1, PSMs_uniprot_df_1, 
					PSMs_extended_df_2, PSMs_ensembl_ref_df_2, PSMs_uniprot_iso_df_2, PSMs_uniprot_df_2, 
					PSMs_extended_df_3, PSMs_ensembl_ref_df_3, PSMs_uniprot_iso_df_3, PSMs_uniprot_df_3, 
					feature_name, annotation_x, plots_directory_path, title_info, level_of_analysis):
	plots = Plots()
	figure_name = plots_directory_path + "/experiments_numOfHits_" + feature_name + title_info + ".png"

	exp1_extended = find_targets_decoys_on_threshold(PSMs_extended_df_1, feature_name, 1, 2)
	exp1_ensembl_ref = find_targets_decoys_on_threshold(PSMs_ensembl_ref_df_1, feature_name, 1, 2)
	exp1_uniprot_iso = find_targets_decoys_on_threshold(PSMs_uniprot_iso_df_1, feature_name, 1, 2)
	exp1_uniprot = find_targets_decoys_on_threshold(PSMs_uniprot_df_1, feature_name, 1, 2)
	
	exp2_extended = find_targets_decoys_on_threshold(PSMs_extended_df_2, feature_name, 1, 2)
	exp2_ensembl_ref = find_targets_decoys_on_threshold(PSMs_ensembl_ref_df_2, feature_name, 1, 2)
	exp2_uniprot_iso = find_targets_decoys_on_threshold(PSMs_uniprot_iso_df_2, feature_name, 1, 2)
	exp2_uniprot = find_targets_decoys_on_threshold(PSMs_uniprot_df_2, feature_name, 1, 2)
	
	exp3_extended = find_targets_decoys_on_threshold(PSMs_extended_df_3, feature_name, 1, 2)
	exp3_ensembl_ref = find_targets_decoys_on_threshold(PSMs_ensembl_ref_df_3, feature_name, 1, 2)
	exp3_uniprot_iso = find_targets_decoys_on_threshold(PSMs_uniprot_iso_df_3, feature_name, 1, 2)
	exp3_uniprot = find_targets_decoys_on_threshold(PSMs_uniprot_df_3, feature_name, 1, 2)

	fig, ax = plt.subplots(figsize=(15,10), dpi= 200)
	matplotlib.rcParams['legend.fontsize'] = 24
	
	xlabel = 'Estimated ' + level_of_analysis + '-level FDR'
	ylabel = 'Number of ' + level_of_analysis + 's'
	
	extendedDB_patch = mpatches.Patch(color='blue', label='Extended DB')
	ensembl_ref_patch = mpatches.Patch(color='green', label='Ensembl-reference')
	uniprot_iso_patch = mpatches.Patch(color='brown', label='UniProt-isoforms')
	uniprotDB_patch = mpatches.Patch(color='orange', label='UniProt DB')
	
	normalizedData = exp1_extended[0]
	
	plot1, = plt.plot(exp1_extended[1], normalizedData, c = 'blue', linewidth=3, label='Experiment P010694')
	
	plot2, = plt.plot(exp1_ensembl_ref[1], exp1_ensembl_ref[0], c = 'green', linewidth=3)
	
	plot3, = plt.plot(exp1_uniprot_iso[1], exp1_uniprot_iso[0], c = 'brown', linewidth=3)
	
	normalizedData = exp1_uniprot[0]
	
	plot4, = plt.plot(exp1_uniprot[1], normalizedData, c = 'orange', linewidth=3)
	
	maxy = max([exp1_extended[0][-1], exp1_ensembl_ref[0][-1], exp1_uniprot_iso[0][-1], exp1_uniprot[0][-1]])*1.02
	total_max = maxy
	plt.annotate('Experiment P010694', xy=(0.042, maxy), fontsize=18)
	
	normalizedData = exp2_extended[0]
	
	plot5, = plt.plot(exp2_extended[1], normalizedData, c = 'blue', linewidth=3, linestyle='--', label='Experiment P010747')
	
	plot6, = plt.plot(exp2_ensembl_ref[1], exp2_ensembl_ref[0], c = 'green', linewidth=3, linestyle='--')
	
	plot7, = plt.plot(exp2_uniprot_iso[1], exp2_uniprot_iso[0], c = 'brown', linewidth=3, linestyle='--')
	
	normalizedData = exp2_uniprot[0]
	
	plot8, = plt.plot(exp2_uniprot[1], normalizedData, c = 'orange', linewidth=3, linestyle='--')
	
	maxy = max([exp2_extended[0][-1], exp2_ensembl_ref[0][-1], exp2_uniprot_iso[0][-1], exp2_uniprot[0][-1]])*1.02
	total_max = max([total_max, maxy])
	plt.annotate('Experiment P010747', xy=(0.042, maxy), fontsize=18)
	
	normalizedData = exp3_extended[0]
	
	plot9, = plt.plot(exp3_extended[1], normalizedData, c = 'blue', linewidth=3, linestyle=':',label='Experiment P013107')
	
	plot10, = plt.plot(exp3_ensembl_ref[1], exp3_ensembl_ref[0], c = 'green', linewidth=3, linestyle=':')
	
	plot11, = plt.plot(exp3_uniprot_iso[1], exp3_uniprot_iso[0], c = 'brown', linewidth=3, linestyle=':')
	
	normalizedData = exp3_uniprot[0]
	
	plot12, = plt.plot(exp3_uniprot[1], normalizedData, c = 'orange', linewidth=3, linestyle=':')
	
	maxy = max([exp3_extended[0][-1], exp3_ensembl_ref[0][-1], exp3_uniprot_iso[0][-1], exp3_uniprot[0][-1]])*1.02
	total_max = max([total_max, maxy])
	plt.annotate('Experiment P013107', xy=(0.042, maxy), fontsize=18)
	
	ax.legend(handles=[extendedDB_patch, ensembl_ref_patch, uniprot_iso_patch, uniprotDB_patch])
	
	plt.ylim(top=total_max+(total_max*0.08))
	plt.xticks(fontsize=27)
	plt.yticks(fontsize=27)
	plt.xlabel(xlabel, fontsize=30, labelpad=20)
	plt.ylabel(ylabel, fontsize=30, labelpad=20)
	
	
	plt.tight_layout()
	plt.savefig(figure_name)
	plt.close()
	
	figure_name = plots_directory_path + "/experiments_numOfHits_" + feature_name + title_info + "_fullYaxis.png"
	
	fig, ax = plt.subplots(figsize=(15,10), dpi= 200)
	matplotlib.rcParams['legend.fontsize'] = 24
	
	normalizedData = exp1_extended[0]
	exp1_extended[1] = [0.0] + exp1_extended[1]
	normalizedData = [0.0] + normalizedData
	plot1, = plt.plot(exp1_extended[1], normalizedData, c = 'blue', linewidth=3, label='Experiment P010694')
	
	exp1_ensembl_ref[1] = [0.0] + exp1_ensembl_ref[1]
	exp1_ensembl_ref[0] = [0.0] + exp1_ensembl_ref[0]
	plot2, = plt.plot(exp1_ensembl_ref[1], exp1_ensembl_ref[0], c = 'green', linewidth=3)
	
	exp1_uniprot_iso[1] = [0.0] + exp1_uniprot_iso[1]
	exp1_uniprot_iso[0] = [0.0] + exp1_uniprot_iso[0]
	plot3, = plt.plot(exp1_uniprot_iso[1], exp1_uniprot_iso[0], c = 'brown', linewidth=3)
	
	normalizedData = exp1_uniprot[0]
	exp1_uniprot[1] = [0.0] + exp1_uniprot[1]
	normalizedData = [0.0] + normalizedData
	plot4, = plt.plot(exp1_uniprot[1], normalizedData, c = 'orange', linewidth=3)
	
	maxy = max([exp1_extended[0][-1], exp1_ensembl_ref[0][-1], exp1_uniprot_iso[0][-1], exp1_uniprot[0][-1]])*1.02
	total_max = maxy
	plt.annotate('Experiment P010694', xy=(0.042, maxy), fontsize=18)
	
	normalizedData = exp2_extended[0]
	exp2_extended[1] = [0.0] + exp2_extended[1]
	normalizedData = [0.0] + normalizedData
	plot5, = plt.plot(exp2_extended[1], normalizedData, c = 'blue', linewidth=3, linestyle='--', label='Experiment P010747')
	
	exp2_ensembl_ref[1] = [0.0] + exp2_ensembl_ref[1]
	exp2_ensembl_ref[0] = [0.0] + exp2_ensembl_ref[0]
	plot6, = plt.plot(exp2_ensembl_ref[1], exp2_ensembl_ref[0], c = 'green', linewidth=3, linestyle='--')
	
	exp2_uniprot_iso[1] = [0.0] + exp2_uniprot_iso[1]
	exp2_uniprot_iso[0] = [0.0] + exp2_uniprot_iso[0]
	plot7, = plt.plot(exp2_uniprot_iso[1], exp2_uniprot_iso[0], c = 'brown', linewidth=3, linestyle='--')
	
	normalizedData = exp2_uniprot[0]
	exp2_uniprot[1] = [0.0] + exp2_uniprot[1]
	normalizedData = [0.0] + normalizedData
	plot8, = plt.plot(exp2_uniprot[1], normalizedData, c = 'orange', linewidth=3, linestyle='--')
	
	maxy = max([exp2_extended[0][-1], exp2_ensembl_ref[0][-1], exp2_uniprot_iso[0][-1], exp2_uniprot[0][-1]])*1.02
	total_max = max([total_max, maxy])
	plt.annotate('Experiment P010747', xy=(0.042, maxy), fontsize=18)
	
	normalizedData = exp3_extended[0]
	exp3_extended[1] = [0.0] + exp3_extended[1]
	normalizedData = [0.0] + normalizedData
	plot9, = plt.plot(exp3_extended[1], normalizedData, c = 'blue', linewidth=3, linestyle=':',label='Experiment P013107')
	
	exp3_ensembl_ref[1] = [0.0] + exp3_ensembl_ref[1]
	exp3_ensembl_ref[0] = [0.0] + exp3_ensembl_ref[0]
	plot10, = plt.plot(exp3_ensembl_ref[1], exp3_ensembl_ref[0], c = 'green', linewidth=3, linestyle=':')
	
	exp3_uniprot_iso[1] = [0.0] + exp3_uniprot_iso[1]
	exp3_uniprot_iso[0] = [0.0] + exp3_uniprot_iso[0]
	plot11, = plt.plot(exp3_uniprot_iso[1], exp3_uniprot_iso[0], c = 'brown', linewidth=3, linestyle=':')
	
	normalizedData = exp3_uniprot[0]
	exp3_uniprot[1] = [0.0] + exp3_uniprot[1]
	normalizedData = [0.0] + normalizedData
	plot12, = plt.plot(exp3_uniprot[1], normalizedData, c = 'orange', linewidth=3, linestyle=':')
	
	maxy = max([exp3_extended[0][-1], exp3_ensembl_ref[0][-1], exp3_uniprot_iso[0][-1], exp3_uniprot[0][-1]])*1.02
	total_max = max([total_max, maxy])
	plt.annotate('Experiment P013107', xy=(0.042, maxy), fontsize=18)
	
	ax.legend(handles=[extendedDB_patch, ensembl_ref_patch, uniprot_iso_patch, uniprotDB_patch])
	
	plt.ylim(0.0, total_max+(total_max*0.08))
	plt.xticks(fontsize=27)
	plt.yticks(fontsize=27)
	plt.xlabel(xlabel, fontsize=30, labelpad=20)
	plt.ylabel(ylabel, fontsize=30, labelpad=20)
	plt.tight_layout()
	plt.savefig(figure_name)
	plt.close()
	
	
def evalue_number_of_hits_Percolator(PSMs_extended_df_1, PSMs_standard_df_1, PSMs_extended_df_2, PSMs_standard_df_2, PSMs_extended_df_3, PSMs_standard_df_3, annotation_x, plots_directory_path, level_of_analysis):
	plots = Plots()
	figure_name = plots_directory_path + "/experiments_numOfHits_Percolator_extendedDB.png"

	exp1_extended = find_targets_decoys_on_threshold(PSMs_extended_df_1, 'perc_qvalue', 1, 4)
	exp1_standard = find_targets_decoys_on_threshold(PSMs_standard_df_1, 'perc_qvalue', 1, 4)
	
	exp2_extended = find_targets_decoys_on_threshold(PSMs_extended_df_2, 'perc_qvalue', 1, 4)
	exp2_standard = find_targets_decoys_on_threshold(PSMs_standard_df_2, 'perc_qvalue', 1, 4)
	
	exp3_extended = find_targets_decoys_on_threshold(PSMs_extended_df_3, 'perc_qvalue', 1, 4)
	exp3_standard = find_targets_decoys_on_threshold(PSMs_standard_df_3, 'perc_qvalue', 1, 4)

	fig, ax = plt.subplots(figsize=(15,10), dpi= 200)
	matplotlib.rcParams['legend.fontsize'] = 24
	
	xlabel = 'Estimated ' + level_of_analysis + '-level FDR'
	ylabel = 'Number of ' + level_of_analysis + 's'
	
	extended_patch = mpatches.Patch(color='mediumblue', label='Percolator extended features')
	standard_patch = mpatches.Patch(color='firebrick', label='Percolator standard features')
	
	plt.plot(exp1_extended[1], exp1_extended[0], c = 'mediumblue', linewidth=4)
	plt.plot(exp1_standard[1], exp1_standard[0], c = 'firebrick', linewidth=4)

	plt.annotate('Experiment P010694', xy=(0.04, exp1_extended[0][-1]*1.02), fontsize=18)
	
	plt.plot(exp2_extended[1], exp2_extended[0], c = 'mediumblue', linewidth=4, linestyle='--')
	plt.plot(exp2_standard[1], exp2_standard[0], c = 'firebrick', linewidth=4, linestyle='--')

	plt.annotate('Experiment P010747', xy=(0.04, exp2_extended[0][-1]*1.03), fontsize=18)
	
	plt.plot(exp3_extended[1], exp3_extended[0], c = 'mediumblue', linewidth=4, linestyle=':')
	plt.plot(exp3_standard[1], exp3_standard[0], c = 'firebrick', linewidth=4, linestyle=':')

	plt.annotate('Experiment P013107', xy=(0.04, exp3_extended[0][-1]*1.04), fontsize=18)
	
	
	miny = min(exp1_standard[0] + exp2_standard[0] + exp3_standard[0])
	maxy = max(exp2_extended[0] + exp1_extended[0] + exp3_extended[0])
	y_range = maxy - miny
	
	plt.ylim(miny-(y_range*0.08),maxy+(y_range*0.08))
	
	plt.xticks(fontsize=27)
	plt.yticks(fontsize=27)
	plt.xlabel(xlabel, fontsize=30, labelpad=20)
	plt.ylabel(ylabel, fontsize=30, labelpad=20)
	plt.legend(handles=[extended_patch, standard_patch])
	plt.tight_layout()
	plt.savefig(figure_name)
	plt.close()
	
	figure_name = plots_directory_path + "/experiments_numOfHits_Percolator_extendedDB_fullYaxis.png"
	
	fig, ax = plt.subplots(figsize=(15,10), dpi= 200)
	matplotlib.rcParams['legend.fontsize'] = 24
	
	plt.plot(exp1_extended[1], exp1_extended[0], c = 'mediumblue', linewidth=4)
	plt.plot(exp1_standard[1], exp1_standard[0], c = 'firebrick', linewidth=4)
	plt.annotate('Experiment P010694', xy=(0.04, exp1_extended[0][-1]*1.02), fontsize=18)
	
	plt.plot(exp2_extended[1], exp2_extended[0], c = 'mediumblue', linewidth=4, linestyle='--')
	plt.plot(exp2_standard[1], exp2_standard[0], c = 'firebrick', linewidth=4, linestyle='--')
	plt.annotate('Experiment P010747', xy=(0.04, exp2_extended[0][-1]*1.03), fontsize=18)
	
	plt.plot(exp3_extended[1], exp3_extended[0], c = 'mediumblue', linewidth=4, linestyle=':')
	plt.plot(exp3_standard[1], exp3_standard[0], c = 'firebrick', linewidth=4, linestyle=':')
	plt.annotate('Experiment P013107', xy=(0.04, exp3_extended[0][-1]*1.04), fontsize=18)
	
	
	miny = min(exp1_standard[0] + exp2_standard[0] + exp3_standard[0])
	maxy = max(exp2_extended[0] + exp1_extended[0] + exp3_extended[0])
	y_range = maxy - miny
	
	plt.ylim(0.0,maxy+(y_range*0.08))
	
	plt.xticks(fontsize=27)
	plt.yticks(fontsize=27)
	plt.xlabel(xlabel, fontsize=30, labelpad=20)
	plt.ylabel(ylabel, fontsize=30, labelpad=20)
	plt.legend(handles=[extended_patch, standard_patch])
	plt.tight_layout()
	plt.savefig(figure_name)
	plt.close()

	
	

P010694_extendedDB_bench_df = pd.read_csv(sys.argv[1])
if 'benchmark_aggreg_q-value' in P010694_extendedDB_bench_df:
	P010694_extendedDB_bench_df = P010694_extendedDB_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P010694_extendedDB_bench_df:
	P010694_extendedDB_bench_df = P010694_extendedDB_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P010694_extendedDB_extra_df = pd.read_csv(sys.argv[2])
if 'extra_features_aggreg_q-value' in P010694_extendedDB_extra_df:
	P010694_extendedDB_extra_df = P010694_extendedDB_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P010694_extendedDB_extra_df:
	P010694_extendedDB_extra_df = P010694_extendedDB_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})


#------------------------------------------

P010747_extendedDB_bench_df = pd.read_csv(sys.argv[3])
P010747_extendedDB_bench_df = P010747_extendedDB_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
if 'benchmark_aggreg_q-value' in P010747_extendedDB_bench_df:
	P010747_extendedDB_bench_df = P010747_extendedDB_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P010747_extendedDB_bench_df:
	P010747_extendedDB_bench_df = P010747_extendedDB_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P010747_extendedDB_extra_df = pd.read_csv(sys.argv[4])
if 'extra_features_aggreg_q-value' in P010747_extendedDB_extra_df:
	P010747_extendedDB_extra_df = P010747_extendedDB_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P010747_extendedDB_extra_df:
	P010747_extendedDB_extra_df = P010747_extendedDB_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})

#------------------------------------------

P013107_extendedDB_bench_df = pd.read_csv(sys.argv[5])
if 'benchmark_aggreg_q-value' in P013107_extendedDB_bench_df:
	P013107_extendedDB_bench_df = P013107_extendedDB_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P013107_extendedDB_bench_df:
	P013107_extendedDB_bench_df = P013107_extendedDB_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P013107_extendedDB_extra_df = pd.read_csv(sys.argv[6])
if 'extra_features_aggreg_q-value' in P013107_extendedDB_extra_df:
	P013107_extendedDB_extra_df = P013107_extendedDB_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P013107_extendedDB_extra_df:
	P013107_extendedDB_extra_df = P013107_extendedDB_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})

#------------------------------------------

P010694_uniprot_bench_df = pd.read_csv(sys.argv[7])
if 'benchmark_aggreg_q-value' in P010694_uniprot_bench_df:
	P010694_uniprot_bench_df = P010694_uniprot_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P010694_uniprot_bench_df:
	P010694_uniprot_bench_df = P010694_uniprot_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010694_uniprot_bench_df['psm_type'][P010694_uniprot_bench_df.Decoy == 1] = "canonical"

P010694_uniprot_extra_df = pd.read_csv(sys.argv[8])
if 'extra_features_aggreg_q-value' in P010694_uniprot_extra_df:
	P010694_uniprot_extra_df = P010694_uniprot_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P010694_uniprot_extra_df:
	P010694_uniprot_extra_df = P010694_uniprot_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010694_uniprot_extra_df['psm_type'][P010694_uniprot_extra_df.Decoy == 1] = "canonical"

#------------------------------------------

P010747_uniprot_bench_df = pd.read_csv(sys.argv[9])
if 'benchmark_aggreg_q-value' in P010747_uniprot_bench_df:
	P010747_uniprot_bench_df = P010747_uniprot_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P010747_uniprot_bench_df:
	P010747_uniprot_bench_df = P010747_uniprot_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010747_uniprot_bench_df['psm_type'][P010747_uniprot_bench_df.Decoy == 1] = "canonical"

P010747_uniprot_extra_df = pd.read_csv(sys.argv[10])
if 'extra_features_aggreg_q-value' in P010747_uniprot_extra_df:
	P010747_uniprot_extra_df = P010747_uniprot_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P010747_uniprot_extra_df:
	P010747_uniprot_extra_df = P010747_uniprot_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010747_uniprot_extra_df['psm_type'][P010747_uniprot_extra_df.Decoy == 1] = "canonical"

#------------------------------------------

P013107_uniprot_bench_df = pd.read_csv(sys.argv[11])
if 'benchmark_aggreg_q-value' in P013107_uniprot_bench_df:
	P013107_uniprot_bench_df = P013107_uniprot_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P013107_uniprot_bench_df:
	P013107_uniprot_bench_df = P013107_uniprot_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P013107_uniprot_bench_df['psm_type'][P013107_uniprot_bench_df.Decoy == 1] = "canonical"

P013107_uniprot_extra_df = pd.read_csv(sys.argv[12])
if 'extra_features_aggreg_q-value' in P013107_uniprot_extra_df:
	P013107_uniprot_extra_df = P013107_uniprot_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P013107_uniprot_extra_df:
	P013107_uniprot_extra_df = P013107_uniprot_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P013107_uniprot_extra_df['psm_type'][P013107_uniprot_extra_df.Decoy == 1] = "canonical"

#------------------------------------------

P010694_ensembl_ref_bench_df = pd.read_csv(sys.argv[13])
if 'benchmark_aggreg_q-value' in P010694_ensembl_ref_bench_df:
	P010694_ensembl_ref_bench_df = P010694_ensembl_ref_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P010694_ensembl_ref_bench_df:
	P010694_ensembl_ref_bench_df = P010694_ensembl_ref_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P010694_ensembl_ref_extra_df = pd.read_csv(sys.argv[14])
if 'extra_features_aggreg_q-value' in P010694_ensembl_ref_extra_df:
	P010694_ensembl_ref_extra_df = P010694_ensembl_ref_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P010694_ensembl_ref_extra_df:
	P010694_ensembl_ref_extra_df = P010694_ensembl_ref_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})

#------------------------------------------

P010747_ensembl_ref_bench_df = pd.read_csv(sys.argv[15])
if 'benchmark_aggreg_q-value' in P010747_ensembl_ref_bench_df:
	P010747_ensembl_ref_bench_df = P010747_ensembl_ref_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P010747_ensembl_ref_bench_df:
	P010747_ensembl_ref_bench_df = P010747_ensembl_ref_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P010747_ensembl_ref_extra_df = pd.read_csv(sys.argv[16])
if 'extra_features_aggreg_q-value' in P010747_ensembl_ref_extra_df:
	P010747_ensembl_ref_extra_df = P010747_ensembl_ref_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P010747_ensembl_ref_extra_df:
	P010747_ensembl_ref_extra_df = P010747_ensembl_ref_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})

#------------------------------------------

P013107_ensembl_ref_bench_df = pd.read_csv(sys.argv[17])
if 'benchmark_aggreg_q-value' in P013107_ensembl_ref_bench_df:
	P013107_ensembl_ref_bench_df = P013107_ensembl_ref_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
elif 'benchmark_q-value' in P013107_ensembl_ref_bench_df:
	P013107_ensembl_ref_bench_df = P013107_ensembl_ref_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep"})

P013107_ensembl_ref_extra_df = pd.read_csv(sys.argv[18])
if 'extra_features_aggreg_q-value' in P013107_ensembl_ref_extra_df:
	P013107_ensembl_ref_extra_df = P013107_ensembl_ref_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
elif 'extra_features_q-value' in P013107_ensembl_ref_extra_df:
	P013107_ensembl_ref_extra_df = P013107_ensembl_ref_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep"})

#------------------------------------------

P010694_uniprot_iso_bench_df = pd.read_csv(sys.argv[19])
if 'benchmark_aggreg_q-value' in P010694_uniprot_iso_bench_df:
	P010694_uniprot_iso_bench_df = P010694_uniprot_iso_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P010694_uniprot_iso_bench_df:
	P010694_uniprot_iso_bench_df = P010694_uniprot_iso_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010694_uniprot_iso_bench_df['psm_type'][P010694_uniprot_iso_bench_df.Decoy == 1] = "canonical"

P010694_uniprot_iso_extra_df = pd.read_csv(sys.argv[20])
if 'extra_features_aggreg_q-value' in P010694_uniprot_iso_extra_df:
	P010694_uniprot_iso_extra_df = P010694_uniprot_iso_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P010694_uniprot_iso_extra_df:
	P010694_uniprot_iso_extra_df = P010694_uniprot_iso_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010694_uniprot_iso_extra_df['psm_type'][P010694_uniprot_iso_extra_df.Decoy == 1] = "canonical"

#------------------------------------------

P010747_uniprot_iso_bench_df = pd.read_csv(sys.argv[21])
if 'benchmark_aggreg_q-value' in P010747_uniprot_iso_bench_df:
	P010747_uniprot_iso_bench_df = P010747_uniprot_iso_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P010747_uniprot_iso_bench_df:
	P010747_uniprot_iso_bench_df = P010747_uniprot_iso_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010747_uniprot_iso_bench_df['psm_type'][P010747_uniprot_iso_bench_df.Decoy == 1] = "canonical"

P010747_uniprot_iso_extra_df = pd.read_csv(sys.argv[22])
if 'extra_features_aggreg_q-value' in P010747_uniprot_iso_extra_df:
	P010747_uniprot_iso_extra_df = P010747_uniprot_iso_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P010747_uniprot_iso_extra_df:
	P010747_uniprot_iso_extra_df = P010747_uniprot_iso_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P010747_uniprot_iso_extra_df['psm_type'][P010747_uniprot_iso_extra_df.Decoy == 1] = "canonical"

#------------------------------------------

P013107_uniprot_iso_bench_df = pd.read_csv(sys.argv[23])
if 'benchmark_aggreg_q-value' in P013107_uniprot_iso_bench_df:
	P013107_uniprot_iso_bench_df = P013107_uniprot_iso_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'benchmark_q-value' in P013107_uniprot_iso_bench_df:
	P013107_uniprot_iso_bench_df = P013107_uniprot_iso_bench_df.rename(columns={"benchmark_score": "perc_score", "benchmark_q-value": "perc_qvalue", "benchmark_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P013107_uniprot_iso_bench_df['psm_type'][P013107_uniprot_iso_bench_df.Decoy == 1] = "canonical"

P013107_uniprot_iso_extra_df = pd.read_csv(sys.argv[24])
if 'extra_features_aggreg_q-value' in P013107_uniprot_iso_extra_df:
	P013107_uniprot_iso_extra_df = P013107_uniprot_iso_extra_df.rename(columns={"extra_features_aggreg_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
elif 'extra_features_q-value' in P013107_uniprot_iso_extra_df:
	P013107_uniprot_iso_extra_df = P013107_uniprot_iso_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_q-value": "perc_qvalue", "extra_features_posterior_error_prob": "perc_pep",  "peptide_type": "psm_type"})
P013107_uniprot_iso_extra_df['psm_type'][P013107_uniprot_iso_extra_df.Decoy == 1] = "canonical"

plots_directory_path = sys.argv[25]

#Should be 'PSM' or 'peptide'
level_of_analysis = sys.argv[26]
#------------------------------------------


EXPERIMENTS = ['P010694', 'P010747', 'P013107']

experiment_dfs = {'P010694': [P010694_extendedDB_bench_df, P010694_extendedDB_extra_df, P010694_uniprot_bench_df, P010694_uniprot_extra_df,
				P010694_ensembl_ref_bench_df, P010694_ensembl_ref_extra_df, P010694_uniprot_iso_bench_df, P010694_uniprot_iso_extra_df],
		  'P010747': [P010747_extendedDB_bench_df, P010747_extendedDB_extra_df, P010747_uniprot_bench_df, P010747_uniprot_extra_df,
		  		P010747_ensembl_ref_bench_df, P010747_ensembl_ref_extra_df, P010747_uniprot_iso_bench_df, P010747_uniprot_iso_extra_df],
		  'P013107': [P013107_extendedDB_bench_df, P013107_extendedDB_extra_df, P013107_uniprot_bench_df, P013107_uniprot_extra_df,
		  		P013107_ensembl_ref_bench_df, P013107_ensembl_ref_extra_df, P013107_uniprot_iso_bench_df, P013107_uniprot_iso_extra_df]}


evalue_number_of_hits_databases(experiment_dfs['P010694'][0], experiment_dfs['P010694'][4], experiment_dfs['P010694'][6], experiment_dfs['P010694'][2],
			 experiment_dfs['P010747'][0], experiment_dfs['P010747'][4], experiment_dfs['P010747'][6], experiment_dfs['P010747'][2],
			  experiment_dfs['P013107'][0], experiment_dfs['P013107'][4], experiment_dfs['P013107'][6], experiment_dfs['P013107'][2], 'engine_score', -7, plots_directory_path, '', level_of_analysis)


evalue_number_of_hits_Percolator(experiment_dfs['P010694'][1], experiment_dfs['P010694'][0], experiment_dfs['P010747'][1], experiment_dfs['P010747'][0], experiment_dfs['P013107'][1], experiment_dfs['P013107'][0], -10, plots_directory_path, level_of_analysis)


