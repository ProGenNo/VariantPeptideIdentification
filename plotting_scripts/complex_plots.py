from plots import *
from RT import *
from matplotlib_venn import venn2, venn3, venn2_circles

class SearchSpace2D:

	def __init__(self, RT):
		self.plots = Plots()
		self.RT = RT
		
	def venn_PSMs(self, bench_PSMs_df, ext_PSMs_df, database_name, plots_directory_path, psm_type):
		figurename = plots_directory_path + '/2D_plots/confident_'+psm_type+'_venn.png'
		transparent = False
		
		standard_accepted = bench_PSMs_df.loc[(bench_PSMs_df['psm_type'] == psm_type) & (bench_PSMs_df['perc_qvalue'] <= 0.01)]['PSMId'].tolist()
		extra_accepted = ext_PSMs_df.loc[(ext_PSMs_df['psm_type'] == psm_type) & (ext_PSMs_df['perc_qvalue'] <= 0.01)]['PSMId'].tolist()
		
		print(len(standard_accepted), len(extra_accepted))
		
		df_standard = pd.DataFrame(standard_accepted,columns=['standard_features'])
		df_extra = pd.DataFrame(extra_accepted,columns=['added_features'])
		
		A = set(df_standard.standard_features)
		B = set(df_extra.added_features)
		
		print(len(A), len(B))
		
		plt.figure(figsize=(15,10), dpi= 200)
		
		ax = plt.gca() 
		
		v = venn2([A, B], set_labels = ('Standard features', 'Extended features'))
		v.get_patch_by_id('10').set_alpha(0.4)
		v.get_patch_by_id('10').set_facecolor('firebrick')
		v.get_patch_by_id('10').set_edgecolor('darkred')
		v.get_patch_by_id('10').set_linewidth(3.0)

		v.get_patch_by_id('01').set_alpha(0.5)
		v.get_patch_by_id('01').set_facecolor('darkslateblue')
		v.get_patch_by_id('01').set_edgecolor('darkblue')
		v.get_patch_by_id('01').set_linewidth(3.0)
		
		c = venn2_circles([A, B], linestyle='-')
		c[0].set_lw(5.0)
		c[0].set_ls('-')
		c[0].set_color('firebrick')
		c[0].set_alpha(0.4)
		c[1].set_lw(5.0)
		c[1].set_ls('-')
		c[1].set_color('darkslateblue')
		c[1].set_alpha(0.5)
		
		out=v 
		for text in out.set_labels:
			text.set_fontsize(30)
		for text in out.subset_labels:
			text.set_fontsize(25)
		for idx, subset in enumerate(out.subset_labels):
			out.subset_labels[idx].set_visible(True)
		for idx, subset in enumerate(out.set_labels):
			out.set_labels[idx].set_visible(True)
		plt.tight_layout()
		plt.savefig(figurename, transparent=transparent)
		plt.close()
		
	def conf_targets_acceptedGroups_RTresiduals_fragAng(self, bench_PSMs_df, ext_PSMs_df, database_name, plots_directory_path, psm_types, perc_feature_set, spectra_feature, ylabel, bins_N):
		
		bench_PSMs_df = bench_PSMs_df[(bench_PSMs_df['perc_qvalue'] <= 0.01) & (bench_PSMs_df['psm_type'].isin(psm_types))]
		ext_PSMs_df = ext_PSMs_df[(ext_PSMs_df['perc_qvalue'] <= 0.01) & (ext_PSMs_df['psm_type'].isin(psm_types))]
		
		both_df = bench_PSMs_df.merge(ext_PSMs_df, how='inner', on='PSMId', suffixes=('', '_duplicate'))
		dupl_cols = both_df.columns[both_df.columns.str.endswith('_duplicate')]
		both_df.drop(dupl_cols, axis=1, inplace=True)
		
		standard_df = bench_PSMs_df.merge(ext_PSMs_df.drop_duplicates(), on=['PSMId'], how='left', indicator=True, suffixes=('', '_duplicate'))
		dupl_cols = standard_df.columns[standard_df.columns.str.endswith('_duplicate')]
		standard_df.drop(dupl_cols, axis=1, inplace=True)
		standard_df = standard_df[standard_df['_merge'] == 'left_only']
		
		extended_df = ext_PSMs_df.merge(bench_PSMs_df.drop_duplicates(), on=['PSMId'], how='left', indicator=True, suffixes=('', '_duplicate'))
		dupl_cols = extended_df.columns[extended_df.columns.str.endswith('_duplicate')]
		extended_df.drop(dupl_cols, axis=1, inplace=True)
		extended_df = extended_df[extended_df['_merge'] == 'left_only']
		
		print(both_df.shape)
		print(standard_df.shape)
		print(extended_df.shape)
		
		y_label = ylabel
		x_label = 'RT error'
		
		all_standard_with_matched_peaks_df = standard_df[standard_df['matched_peaks'] != 0.0]
		number_of_PSM_without_matched_peaks = len(standard_df.index) - len(all_standard_with_matched_peaks_df.index)
		print('Targets with no matched peaks:', number_of_PSM_without_matched_peaks)
		print('Targets min ang sim:', min(all_standard_with_matched_peaks_df[spectra_feature].tolist()))
			
		figurename = plots_directory_path + '/2D_plots' +'/density_' + 'AcceptedStandard' + '_targets_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name + '_' + str(number_of_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(all_standard_with_matched_peaks_df['RT_residuals'].tolist(),
			all_standard_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		
		all_both_with_matched_peaks_df = both_df[both_df['matched_peaks'] != 0.0]
		number_of_PSM_without_matched_peaks = len(both_df.index) - len(all_both_with_matched_peaks_df.index)
		print('Targets with no matched peaks:', number_of_PSM_without_matched_peaks)
		print('Targets min ang sim:', min(all_both_with_matched_peaks_df[spectra_feature].tolist()))
			
		figurename = plots_directory_path + '/2D_plots' +'/density_' + 'AcceptedBoth' + '_targets_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name + '_' + str(number_of_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(all_both_with_matched_peaks_df['RT_residuals'].tolist(),
			all_both_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		
		all_extended_with_matched_peaks_df = extended_df[extended_df['matched_peaks'] != 0.0]
		number_of_PSM_without_matched_peaks = len(extended_df.index) - len(all_extended_with_matched_peaks_df.index)
		print('Targets with no matched peaks:', number_of_PSM_without_matched_peaks)
		print('Targets min ang sim:', min(all_extended_with_matched_peaks_df[spectra_feature].tolist()))
		
		figurename = plots_directory_path + '/2D_plots' +'/density_' + 'AcceptedExtended' + '_targets_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name + '_' + str(number_of_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(all_extended_with_matched_peaks_df['RT_residuals'].tolist(),
			all_extended_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		f = open(plots_directory_path + '/2D_plots' + '/dummy_file_all_targets_acceptedGroups_heatmaps.txt','a+')
		f.write(' ')
		f.close()
		
			
	def all_TD_heatmaps(self, PSMs_df, database_name, plots_directory_path, perc_feature_set, spectra_feature, ylabel, bins_N):
		transparent = False
		figurename = ''
		data_df = None
			
		y_label = ylabel
		x_label = 'RT error'
		
		all_targets_df = PSMs_df[PSMs_df['Decoy'] == 1]
		all_targets_with_matched_peaks_df = all_targets_df[all_targets_df['matched_peaks'] != 0.0]
		number_of_PSM_without_matched_peaks = len(all_targets_df.index) - len(all_targets_with_matched_peaks_df.index)
		print('Targets with no matched peaks:', number_of_PSM_without_matched_peaks)
		print('Targets min ang sim:', min(all_targets_with_matched_peaks_df[spectra_feature].tolist()))
			
		figurename = plots_directory_path + '/2D_plots' +'/density_alltargets_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name + '_' + str(number_of_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(all_targets_with_matched_peaks_df['RT_residuals'].tolist(),
			all_targets_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		decoys_df = PSMs_df[PSMs_df['Decoy'] == -1]
		decoys_with_matched_peaks_df = decoys_df[decoys_df['matched_peaks'] != 0.0]
		number_of_decoy_PSM_without_matched_peaks = len(decoys_df.index) - len(decoys_with_matched_peaks_df.index)
		print('Decoys with no matched peaks:', number_of_decoy_PSM_without_matched_peaks)
		print('Decoys min ang sim:', min(decoys_with_matched_peaks_df[spectra_feature].tolist()))
			
		figurename = plots_directory_path + '/2D_plots' +'/density_decoys_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name +  '_' + str(number_of_decoy_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(decoys_with_matched_peaks_df['RT_residuals'].tolist(),
			decoys_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		f = open(plots_directory_path + '/2D_plots' + '/dummy_file_all_TD_heatmaps.txt','a+')
		f.write(' ')
		f.close()
			
	def all_TD_psmType_heatmaps(self, PSMs_df, database_name, plots_directory_path, perc_feature_set, spectra_feature, ylabel, bins_N, psm_type):
		transparent = False
		figurename = ''
		data_df = None
			
		y_label = ylabel
		x_label = 'RT error'
		
		all_psms_df = PSMs_df[(PSMs_df['Decoy'] == 1) & (PSMs_df['psm_type'] == psm_type)]
		all_psms_with_matched_peaks_df = all_psms_df[all_psms_df['matched_peaks'] != 0.0]
		number_of_PSM_without_matched_peaks = len(all_psms_df.index) - len(all_psms_with_matched_peaks_df.index)
		print('Targets with no matched peaks:', number_of_PSM_without_matched_peaks)
		print('Targets min ang sim:', min(all_psms_with_matched_peaks_df[spectra_feature].tolist()))
			
		figurename = plots_directory_path + '/2D_plots' +'/density_' + psm_type + '_' + spectra_feature + '_bins' + str(bins_N) +'_' + database_name + '_' + str(number_of_PSM_without_matched_peaks) + '_NoMatchedPeaks.png'
		bins = (bins_N,bins_N)
		self.plots.simple_2D_densityplot(all_psms_with_matched_peaks_df['RT_residuals'].tolist(),
			all_psms_with_matched_peaks_df[spectra_feature].tolist(), bins, figurename, x_label, y_label)
			
		f = open(plots_directory_path + '/2D_plots' + '/dummy_file_all_psm_types_heatmaps.txt','a+')
		f.write(' ')
		f.close()

	def conf_variants_scatter_RTresiduals_fragAng(self, bench_PSMs_df, ext_PSMs_df, database_name, plots_directory_path, psm_types):
		figurename = plots_directory_path + '/2D_plots' + "/scatter_RTres_fragmAng_confident_" + '_'.join(psm_types) + ".png"
		
		bench_PSMs_df = bench_PSMs_df[(bench_PSMs_df['perc_qvalue'] <= 0.01) & (bench_PSMs_df['psm_type'].isin(psm_types))]
		ext_PSMs_df = ext_PSMs_df[(ext_PSMs_df['perc_qvalue'] <= 0.01) & (ext_PSMs_df['psm_type'].isin(psm_types))]
		
		both_df = bench_PSMs_df.merge(ext_PSMs_df, how='inner', on='PSMId', suffixes=('', '_duplicate'))
		dupl_cols = both_df.columns[both_df.columns.str.endswith('_duplicate')]
		both_df.drop(dupl_cols, axis=1, inplace=True)
		
		standard_df = bench_PSMs_df.merge(ext_PSMs_df.drop_duplicates(), on=['PSMId'], how='left', indicator=True, suffixes=('', '_duplicate'))
		dupl_cols = standard_df.columns[standard_df.columns.str.endswith('_duplicate')]
		standard_df.drop(dupl_cols, axis=1, inplace=True)
		standard_df = standard_df[standard_df['_merge'] == 'left_only']
		
		extended_df = ext_PSMs_df.merge(bench_PSMs_df.drop_duplicates(), on=['PSMId'], how='left', indicator=True, suffixes=('', '_duplicate'))
		dupl_cols = extended_df.columns[extended_df.columns.str.endswith('_duplicate')]
		extended_df.drop(dupl_cols, axis=1, inplace=True)
		extended_df = extended_df[extended_df['_merge'] == 'left_only']
		
		print(both_df.shape)
		print(standard_df.shape)
		print(extended_df.shape)
		
		x1 = standard_df['RT_residuals'].tolist()
		y1 = standard_df['spectra_angular_similarity'].tolist()
		c1 = 'firebrick'
		s1 = 10
		label1 = 'Standard features'
		x2 = both_df['RT_residuals'].tolist()
		y2 = both_df['spectra_angular_similarity'].tolist()
		c2 = 'purple'
		s2 = 10
		label2 = 'Both feature types'
		x3 = extended_df['RT_residuals'].tolist()
		y3 = extended_df['spectra_angular_similarity'].tolist()
		c3 = 'darkslateblue'
		s3 = 10
		label3 = 'Extended features'
		x_label = 'RT error'
		y_label = 'Spectra angular similarity'
		title = ''
		transparent = False
		
		self.plots.triple_scatterplot(figurename, x1, y1, c1, s1, label1, x2, y2, c2, s2, label2, x3, y3, c3, s3, label3, x_label, y_label, title, transparent)
		
		figurename = plots_directory_path + '/2D_plots' + "/scatter_RTres_fragmAng_confident_standard_"+'_'.join(psm_types)+".png"
		self.plots.mono_scatterplot(figurename, x1, y1, c1, s1, label1, x_label, y_label, title, transparent)
		
		figurename = plots_directory_path + '/2D_plots' + "/scatter_RTres_fragmAng_confident_both_"+'_'.join(psm_types)+".png"
		self.plots.mono_scatterplot(figurename, x2, y2, c2, s2, label2, x_label, y_label, title, transparent)
		
		figurename = plots_directory_path + '/2D_plots' + "/scatter_RTres_fragmAng_confident_extended_"+'_'.join(psm_types)+".png"
		self.plots.mono_scatterplot(figurename, x3, y3, c3, s3, label3, x_label, y_label, title, transparent)
