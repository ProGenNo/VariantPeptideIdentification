import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kde
import numpy as np
import matplotlib as mpl

sns.color_palette("viridis", as_cmap=True)

class Plots:

	def __init__(self):
		print()
		
	def error_dist_target_decoys(self, data_list_1, data_list_2, figurename, title, x_label, label1, label2):
		density1 = kde.gaussian_kde(data_list_1)
		density2 = kde.gaussian_kde(data_list_2)
		
		x1 = np.linspace(min(data_list_1),max(data_list_1),500)
		x2 = np.linspace(min(data_list_2),max(data_list_2),500)

		y1 = density1(x1)
		y2 = density2(x2)
		
		plt.figure(figsize=(15,10), dpi= 200)
		
		plt.plot(x1, y1, c='blue', linewidth=2, label=label1)
		plt.plot(x2, y2, c='green', linewidth=2,label=label2)

		mpl.rcParams['legend.fontsize'] = 28
		plt.xticks(fontsize=25)
		plt.yticks(fontsize=25)
		plt.title(title, fontsize=28)
		plt.xlabel(x_label, fontsize=28)
		plt.ylabel('Density', fontsize=28)
		plt.legend()
		plt.tight_layout()
		plt.savefig(figurename)
		plt.close()
		
	def mono_scatterplot(self, figurename, x1, y1, c1, s1, label1, x_label, y_label, title, transparent):
		self.x_min = min(x1)
		self.x_max = max(x1)
		self.y_min = min(y1)
		self.y_max = max(y1)
		
		y_margin = (self.y_max - self.y_min) * 0.02
		x_margin =  (self.x_max - self.x_min) * 0.02
		self.y_min -= y_margin
		self.y_max += y_margin
		self.x_max += x_margin
		self.x_min -= x_margin
	
		plt.figure(figsize=(15,10), dpi= 200)
		
		plt.scatter(x1, y1, c=c1, s=s1, label=label1)
		
		mpl.rcParams['legend.fontsize'] = 28
		plt.ylim(self.y_min, self.y_max)
		plt.xlim(-6100,6100)
		plt.xticks(fontsize=27)
		plt.yticks(fontsize=27)
		plt.xlabel(x_label, fontsize=30, labelpad=20)
		plt.ylabel(y_label, fontsize=30, labelpad=20)
		plt.title(title, fontsize=20)
		plt.legend()
		plt.tight_layout()
		plt.savefig(figurename, transparent=transparent)
		plt.close()
		
	def triple_scatterplot(self, figurename, x1, y1, c1, s1, label1, x2, y2, c2, s2, label2, x3, y3, c3, s3, label3, x_label, y_label, title, transparent):
		self.x_min = min(x1+x2+x3)
		self.x_max = max(x1+x2+x3)
		self.y_min = min(y1+y2+y3)
		self.y_max = max(y1+y2+y3)
		
		y_margin = (self.y_max - self.y_min) * 0.02
		x_margin =  (self.x_max - self.x_min) * 0.02
		self.y_min -= y_margin
		self.y_max += y_margin
		self.x_max += x_margin
		self.x_min -= x_margin
	
		plt.figure(figsize=(15,10), dpi= 200)
		
		plt.scatter(x1, y1, c=c1, s=s1, label=label1)
		plt.scatter(x2, y2, c=c2, s=s2, label=label2)
		plt.scatter(x3, y3, c=c3, s=s3, label=label3)
		plt.ylim(self.y_min, self.y_max)
		plt.xlim(self.x_min,self.x_max)
		plt.xticks(fontsize=27)
		plt.yticks(fontsize=27)
		plt.xlabel(x_label, fontsize=30, labelpad=20)
		plt.ylabel(y_label, fontsize=30, labelpad=20)
		plt.title(title, fontsize=20)
		plt.legend()
		plt.tight_layout()
		plt.savefig(figurename, transparent=transparent)
		plt.close()
		
	def simple_2D_densityplot(self, x, y, bins, figurename, xlabel, ylabel):
		self.x_min = min(x)
		self.x_max = max(x)
		self.y_min = min(y)
		self.y_max = max(y)
		y_margin = (self.y_max - self.y_min) * 0.02
		x_margin =  (self.x_max - self.x_min) * 0.02
		self.y_min -= y_margin
		self.y_max += y_margin
		self.x_max += x_margin
		self.x_min -= x_margin
		
		transparent = False
		y_label = ylabel
		x_label = xlabel
	
		plt.figure(figsize=(15,10), dpi= 200)
		
		plt.hist2d(x, y, bins=bins, cmap='plasma', cmin = 1)
		
		plt.ylim(self.y_min, self.y_max)
		plt.xlim(-7100,7100)
		plt.xticks(fontsize=27)
		plt.yticks(fontsize=27)
		plt.xlabel(x_label, fontsize=30, labelpad=20)
		plt.ylabel(y_label, fontsize=30, labelpad=20)
		cb=plt.colorbar()
		cb.ax.tick_params(labelsize=27) 
		plt.tight_layout()
		plt.savefig(figurename, transparent=transparent)
		plt.close()
		
	def variants_scatterplot(self, figurename, data_df, x_feature_name, y_feature_name, rt_scatter_flag, x_range, percentiles95, percentiles99, linear_model, database_name, x_label, y_label, title, transparent):
	
		self.x_min = min(data_df[x_feature_name].tolist())
		self.x_max = max(data_df[x_feature_name].tolist())
		self.y_min = min(data_df[y_feature_name].tolist())
		self.y_max = max(data_df[y_feature_name].tolist())
		
		y_margin = (self.y_max - self.y_min) * 0.02
		x_margin =  (self.x_max - self.x_min) * 0.02
		self.y_min -= y_margin
		self.y_max += y_margin
		self.x_max += x_margin
		self.x_min -= x_margin

		df_not_variant = data_df[data_df['psm_type'] == 'canonical']
		df_variant = data_df[data_df['psm_type'] == 'variant']
		
		plt.figure(figsize=(15,10), dpi= 200)
		
		sns.scatterplot(data=df_not_variant, x=x_feature_name, y=y_feature_name, hue="psm_type", s=15, palette="crest")
		sns.scatterplot(data=df_variant, x=x_feature_name, y=y_feature_name, hue="psm_type", s=25, palette="magma")
		
		if rt_scatter_flag==1:
			pos_percentile_95_y = linear_model[0]*x_range + linear_model[1] + percentiles95[0]
			pos_percentile_99_y = linear_model[0]*x_range + linear_model[1] + percentiles99[0]
			neg_percentile_95_y = linear_model[0]*x_range + linear_model[1] + percentiles95[1]
			neg_percentile_99_y = linear_model[0]*x_range + linear_model[1] + percentiles99[1]
			linear_model_y = linear_model[0]*x_range + linear_model[1]
			plt.plot(x_range, linear_model_y, '--', c='gray', label='y=' + "{:.2f}".format(linear_model[0]) +'*x + ' +  "{:.2f}".format(linear_model[1]))
			plt.plot(x_range, pos_percentile_95_y, '--', c='red', label='95th percentile')
			plt.plot(x_range, neg_percentile_95_y, '--', c='red')
			plt.plot(x_range, pos_percentile_99_y, '--', c='darkred',  label='99th percentile')
			plt.plot(x_range, neg_percentile_99_y, '--', c='darkred')
		plt.ylim(self.y_min, self.y_max)
		plt.xlim(-7100,7100)
		plt.xlabel(x_label, fontsize=30 ,labelpad=20)
		plt.ylabel(y_label, fontsize=30, labelpad=20)
		plt.xticks(fontsize=27)
		plt.yticks(fontsize=27)
		plt.title(title, fontsize=25)
		plt.legend(loc='upper left', fontsize=27)
		plt.tight_layout()
		plt.savefig(figurename, transparent=transparent)
		plt.close()

	def violin_plot_4categories(self, datalist, figurename, title, xlabel, ylabel, xticklabels, color, annotate_mean):
		data = datalist
		fig, ax = plt.subplots(figsize=(15,10), dpi= 200)
		quants = [[0.05, 0.95] for i in range(len(data))]
		
		parts = ax.violinplot(data, showmeans=True, showmedians=False, quantiles=quants)
		
		for pc in parts['bodies']:
    			pc.set_facecolor(color)
    			pc.set_alpha(1)
    			
		if annotate_mean ==1:
    			for i in range(len(data)):
    				mean_value = np.mean(np.array(data[i]))
    				plt.text((i+1.15), (mean_value-.01), str(round(mean_value, 3)), fontsize = 14)

		ax.set_title(title)
		ax.set_xlabel(xlabel, fontsize=28, labelpad=20)
		ax.set_ylabel(ylabel, fontsize=28, labelpad=20)
		ax.set_xticks([1.0,2.0,3.0,4.0])
		ax.set_xticklabels(xticklabels)
		ax.yaxis.grid(True)
		plt.xticks(fontsize=26)
		plt.yticks(fontsize=26)
		plt.tight_layout()
		plt.savefig(figurename, transparent=False)
		plt.close()
		
	def violin_plot_8categories(self, datalist, figurename, title, xlabel, ylabel, xticklabels, color, annotate_mean):
		data = datalist
		fig, ax = plt.subplots(figsize=(20,15), dpi= 200)
		quants = [[0.05, 0.95] for i in range(len(data))]
		
		parts = ax.violinplot(data, showmeans=True, showmedians=False, quantiles=quants)
		
		for pc in parts['bodies']:
    			pc.set_facecolor(color)
    			pc.set_alpha(1)
    			
		if annotate_mean ==1:
    			for i in range(len(data)):
    				mean_value = np.mean(np.array(data[i]))
    				plt.text((i+1.15), (mean_value-.01), str(round(mean_value, 3)), fontsize = 14)

		ax.set_title(title)
		ax.set_xlabel(xlabel, fontsize=29, labelpad=10)
		ax.set_ylabel(ylabel, fontsize=29, labelpad=10)
		ax.set_xticks([1.0,2.0,3.0,4.0, 5.0, 6.0, 7.0, 8.0])
		ax.set_xticklabels(xticklabels)
		ax.yaxis.grid(True)
		plt.xticks(fontsize=26, rotation = 60)
		plt.yticks(fontsize=26)
		plt.tight_layout()
		plt.savefig(figurename, transparent=False)
		plt.close()
