import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
from scipy.stats import kde
import sys

from RT import *
from plots import *

def QQ_jointplots(PSMs_extended_df, PSMs_uniprot_df, feature_name, label, plots_directory_path, outliers_info):

	matplotlib.rcParams['xtick.labelsize'] = 14 
	matplotlib.rcParams['ytick.labelsize'] = 14 
	
	extended_targets = PSMs_extended_df.loc[(PSMs_extended_df['Decoy'] == 1)][feature_name].tolist()
	uni_targets = PSMs_uniprot_df.loc[(PSMs_uniprot_df['Decoy'] == 1)][feature_name].tolist()
	extended_decoys = PSMs_extended_df.loc[(PSMs_extended_df['Decoy'] == -1)][feature_name].tolist()
	uni_decoys = PSMs_uniprot_df.loc[(PSMs_uniprot_df['Decoy'] == -1)][feature_name].tolist()
	
	ext_targets_quants = list(np.percentile(extended_targets, range(100)))
	ext_decoys_quants = list(np.percentile(extended_decoys, range(100)))
	uni_targets_quants = list(np.percentile(uni_targets, range(100)))
	uni_decoys_quants = list(np.percentile(uni_decoys, range(100)))
	
	
	
	labels = ['Decoy' for i in range(100)] + ['Target' for i in range(100)]
	extended_data =  ext_decoys_quants + ext_targets_quants
	uniprot_data =  uni_decoys_quants + uni_targets_quants

	sorted_uni_for_line = sorted(uniprot_data)
	
	plot_df = pd.DataFrame()
	plot_df['Extended DB'] = extended_data
	plot_df['UniProt DB'] = uniprot_data
	plot_df['PSM type'] = labels
	
	
	figure_name = plots_directory_path + "/" + feature_name + "_quants_jointplots" + outliers_info + ".png"
	
	fig=plt.figure(figsize=(18, 18))
	ax = fig.add_subplot(111)
	g=sns.jointplot(x="UniProt DB", 
		        y="Extended DB", 
		        hue="PSM type",
		        data=plot_df,
		        linewidth=0,
		        s=45,
		        palette='plasma')
	x0, x1 = g.ax_joint.get_xlim()
	y0, y1 = g.ax_joint.get_ylim()
	lims = [max(x0, y0), min(x1, y1)]
	g.ax_joint.plot(lims, lims, color='black', label='Line of perfect fit')
	
	plt.xlabel(label + ', UniProt DB', fontsize=16, labelpad=20)
	plt.ylabel(label + ', Extended DB', fontsize=16, labelpad=20)
	plt.tight_layout()
	plt.savefig(figure_name,
		            format='png',dpi=200)
	plt.close()
	

def evalue_TD_features_violin_plot(PSMs_extended_df, PSMs_uniprot_df, PSMs_UniProt_isoforms_bench_df, PSMs_Ensembl_reference_bench_df, feature_name, y_label, plots_directory_path):
	plots = Plots()
	
	color = 'indianred'
	
	annotate_mean = 0
	
	figure_name = plots_directory_path + "/" + feature_name + "_databases_TD_violin.png"
	
	extended_targets = PSMs_extended_df.loc[(PSMs_extended_df['Decoy'] == 1)][feature_name].tolist()
	uni_targets = PSMs_uniprot_df.loc[(PSMs_uniprot_df['Decoy'] == 1)][feature_name].tolist()
	extended_decoys = PSMs_extended_df.loc[(PSMs_extended_df['Decoy'] == -1)][feature_name].tolist()
	uni_decoys = PSMs_uniprot_df.loc[(PSMs_uniprot_df['Decoy'] == -1)][feature_name].tolist()
	
	ens_ref_targets = PSMs_Ensembl_reference_bench_df.loc[(PSMs_Ensembl_reference_bench_df['Decoy'] == 1)][feature_name].tolist()
	uni_iso_targets = PSMs_UniProt_isoforms_bench_df.loc[(PSMs_UniProt_isoforms_bench_df['Decoy'] == 1)][feature_name].tolist()
	ens_ref_decoys = PSMs_Ensembl_reference_bench_df.loc[(PSMs_Ensembl_reference_bench_df['Decoy'] == -1)][feature_name].tolist()
	uni_iso_decoys = PSMs_UniProt_isoforms_bench_df.loc[(PSMs_UniProt_isoforms_bench_df['Decoy'] == -1)][feature_name].tolist()
		
	datalist = [extended_targets, ens_ref_targets, uni_iso_targets, uni_targets, extended_decoys, ens_ref_decoys, uni_iso_decoys, uni_decoys]

	xticklabels = ['Targets Extended DB', 'Targets Ensembl-reference', 'Targets UniProt-isoforms', 'Targets UniProtDB', 'Decoys Extended DB', 'Decoys Ensembl-reference', 'Decoys UniProt-isoforms', 'Decoys UniprotDB']
		
	plots.violin_plot_8categories(datalist, figure_name, '', 'PSMs', y_label, xticklabels, color, annotate_mean)
	
	

PSMs_extendedDB_bench_df = pd.read_csv(sys.argv[1])
PSMs_extendedDB_bench_df = PSMs_extendedDB_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})

PSMs_extendedDB_extra_df = pd.read_csv(sys.argv[2])
PSMs_extendedDB_extra_df = PSMs_extendedDB_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
#---------------------------------------------------------

PSMs_uniprot_bench_df = pd.read_csv(sys.argv[3])
PSMs_uniprot_bench_df = PSMs_uniprot_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep", "peptide_type": "psm_type"})
PSMs_uniprot_bench_df['psm_type'][PSMs_uniprot_bench_df.Decoy == 1] = "canonical"

PSMs_uniprot_extra_df = pd.read_csv(sys.argv[4])
PSMs_uniprot_extra_df = PSMs_uniprot_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep", "peptide_type": "psm_type"})
PSMs_uniprot_extra_df['psm_type'][PSMs_uniprot_extra_df.Decoy == 1] = "canonical"
#-----------------------------------------------------------

PSMs_UniProt_isoforms_bench_df = pd.read_csv(sys.argv[5])
PSMs_UniProt_isoforms_bench_df = PSMs_UniProt_isoforms_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep", "peptide_type": "psm_type"})
PSMs_UniProt_isoforms_bench_df['psm_type'][PSMs_UniProt_isoforms_bench_df.Decoy == 1] = "canonical"

PSMs_UniProt_isoforms_extra_df = pd.read_csv(sys.argv[6])
PSMs_UniProt_isoforms_extra_df = PSMs_UniProt_isoforms_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep", "peptide_type": "psm_type"})
PSMs_UniProt_isoforms_extra_df['psm_type'][PSMs_UniProt_isoforms_extra_df.Decoy == 1] = "canonical"
#-----------------------------------------------------------

PSMs_Ensembl_reference_bench_df = pd.read_csv(sys.argv[7])
PSMs_Ensembl_reference_bench_df = PSMs_Ensembl_reference_bench_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})
PSMs_Ensembl_reference_bench_df['psm_type'][PSMs_Ensembl_reference_bench_df.Decoy == 1] = "canonical"

PSMs_Ensembl_reference_extra_df = pd.read_csv(sys.argv[8])
PSMs_Ensembl_reference_extra_df = PSMs_Ensembl_reference_extra_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})
PSMs_Ensembl_reference_extra_df['psm_type'][PSMs_Ensembl_reference_extra_df.Decoy == 1] = "canonical"
#-----------------------------------------------------------

plots_directory_path = sys.argv[9]

rt = RT()
PSMs_extendedDB_bench_df, PSMs_extendedDB_extra_df = rt.find_linear_model(PSMs_extendedDB_bench_df, PSMs_extendedDB_extra_df)

PSMs_extendedDB_bench_df['RT_residuals'].mask(PSMs_extendedDB_bench_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_extendedDB_bench_df['RT_residuals'].mask(PSMs_extendedDB_bench_df['RT_residuals'] > 7000, 7000, inplace=True)
PSMs_extendedDB_extra_df['RT_residuals'].mask(PSMs_extendedDB_extra_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_extendedDB_extra_df['RT_residuals'].mask(PSMs_extendedDB_extra_df['RT_residuals'] > 7000, 7000, inplace=True)

rt = RT()
PSMs_uniprot_bench_df, PSMs_uniprot_extra_df = rt.find_linear_model(PSMs_uniprot_bench_df, PSMs_uniprot_extra_df)
PSMs_uniprot_bench_df['RT_residuals'].mask(PSMs_uniprot_bench_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_uniprot_bench_df['RT_residuals'].mask(PSMs_uniprot_bench_df['RT_residuals'] > 7000, 7000, inplace=True)
PSMs_uniprot_extra_df['RT_residuals'].mask(PSMs_uniprot_extra_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_uniprot_extra_df['RT_residuals'].mask(PSMs_uniprot_extra_df['RT_residuals'] > 7000, 7000, inplace=True)

rt = RT()
PSMs_UniProt_isoforms_bench_df, PSMs_UniProt_isoforms_extra_df = rt.find_linear_model(PSMs_UniProt_isoforms_bench_df, PSMs_UniProt_isoforms_extra_df)
PSMs_UniProt_isoforms_bench_df['RT_residuals'].mask(PSMs_UniProt_isoforms_bench_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_UniProt_isoforms_bench_df['RT_residuals'].mask(PSMs_UniProt_isoforms_bench_df['RT_residuals'] > 7000, 7000, inplace=True)
PSMs_UniProt_isoforms_extra_df['RT_residuals'].mask(PSMs_UniProt_isoforms_extra_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_UniProt_isoforms_extra_df['RT_residuals'].mask(PSMs_UniProt_isoforms_extra_df['RT_residuals'] > 7000, 7000, inplace=True)

rt = RT()
PSMs_Ensembl_reference_bench_df, PSMs_Ensembl_reference_extra_df = rt.find_linear_model(PSMs_Ensembl_reference_bench_df, PSMs_Ensembl_reference_extra_df)
PSMs_Ensembl_reference_bench_df['RT_residuals'].mask(PSMs_Ensembl_reference_bench_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_Ensembl_reference_bench_df['RT_residuals'].mask(PSMs_Ensembl_reference_bench_df['RT_residuals'] > 7000, 7000, inplace=True)
PSMs_Ensembl_reference_extra_df['RT_residuals'].mask(PSMs_Ensembl_reference_extra_df['RT_residuals'] < -7000, -7000, inplace=True)
PSMs_Ensembl_reference_extra_df['RT_residuals'].mask(PSMs_Ensembl_reference_extra_df['RT_residuals'] > 7000, 7000, inplace=True)

#-----------------------------------------------------------

evalue_TD_features_violin_plot(PSMs_extendedDB_bench_df, PSMs_uniprot_bench_df, PSMs_UniProt_isoforms_bench_df, PSMs_Ensembl_reference_bench_df, 'RT_residuals', 'Retention time error', plots_directory_path)

evalue_TD_features_violin_plot(PSMs_extendedDB_bench_df, PSMs_uniprot_bench_df, PSMs_UniProt_isoforms_bench_df, PSMs_Ensembl_reference_bench_df, 'spectra_angular_similarity', 'Spectra angular similarity', plots_directory_path)


##### Jointplots of quantiles targets vs decoys #####

outliers_info = ''

feature_name = 'RT_residuals'
label = 'Quantiles of RT error'
QQ_jointplots(PSMs_extendedDB_bench_df, PSMs_uniprot_bench_df, feature_name, label, plots_directory_path, outliers_info)

feature_name = 'spectra_angular_similarity'
label = 'Quantiles of spectra angular similarity'
QQ_jointplots(PSMs_extendedDB_bench_df, PSMs_uniprot_bench_df, feature_name, label, plots_directory_path, outliers_info)

confident_PSMs_extendedDB_bench_df = PSMs_extendedDB_bench_df[PSMs_extendedDB_bench_df['perc_qvalue'] <= 0.05]

confident_PSMs_uniprot_bench_df = PSMs_uniprot_bench_df[PSMs_uniprot_bench_df['perc_qvalue'] <= 0.05]

outliers_info = '_5FDR'
feature_name = 'engine_score'
label = 'Quantiles of search engine score'
QQ_jointplots(confident_PSMs_extendedDB_bench_df, confident_PSMs_uniprot_bench_df, feature_name, label, plots_directory_path ,outliers_info)



