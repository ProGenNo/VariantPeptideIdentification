import sys
import pandas as pd
from matplotlib_venn import venn2, venn3, venn2_circles
import matplotlib.pyplot as plt
from plots import *

class PeptidePlots:

	def __init__(self):
		self.plots = Plots()
		
	def venn_peptides(self, bench_peptides_df, extra_peptides_df, database_name, plots_directory_path, psm_type):
		figurename = plots_directory_path + '/Peptide_plots/confident_'+psm_type+'_peptides_venn.png'
		transparent = False
		
		standard_accepted = bench_peptides_df.loc[(bench_peptides_df['psm_type'] == psm_type) & (bench_peptides_df['perc_qvalue'] <= 0.01)]['SequenceWithMods'].tolist()
		extra_accepted = extra_peptides_df.loc[(extra_peptides_df['psm_type'] == psm_type) & (extra_peptides_df['perc_qvalue'] <= 0.01)]['SequenceWithMods'].tolist()
		
		print(len(standard_accepted), len(extra_accepted))
		
		df_standard = pd.DataFrame(standard_accepted,columns=['standard_features'])
		df_extra = pd.DataFrame(extra_accepted,columns=['added_features'])
		
		A = set(df_standard.standard_features)
		B = set(df_extra.added_features)
		
		print(len(A), len(B))
		
		plt.figure(figsize=(15,10), dpi= 200)
		
		ax = plt.gca() 
		#out=venn2([A, B], set_labels = ('Standard features', 'Extended features'))
		
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
