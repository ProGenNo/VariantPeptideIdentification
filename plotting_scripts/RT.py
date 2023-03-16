from plots import *
from sklearn.linear_model import LinearRegression
import numpy as np
import pandas as pd

class RT:

	def __init__(self):
		self.plots = Plots()
		
	def find_linear_model(self, bench_PSMs_df, ext_PSMs_df):
		###Linear model####
		
		combined_PSMs_df = pd.concat([bench_PSMs_df, ext_PSMs_df])
		
		#y = combined_PSMs_df.loc[(combined_PSMs_df['psm_type'] == 'canonical') & (combined_PSMs_df['perc_qvalue'] <= 0.01)]['measured_rt'].tolist()
		#X = np.array([combined_PSMs_df.loc[(combined_PSMs_df['psm_type'] == 'canonical') & (combined_PSMs_df['perc_qvalue'] <= 0.01)]['predicted_rt'].tolist()])
		#X = np.transpose(X)
		#print(X.shape)
		
		y = bench_PSMs_df.loc[(bench_PSMs_df['Decoy'] == 1) & (bench_PSMs_df['perc_qvalue'] <= 0.01)]['measured_rt'].tolist()
		X = np.array([bench_PSMs_df.loc[(bench_PSMs_df['Decoy'] == 1) & (bench_PSMs_df['perc_qvalue'] <= 0.01)]['predicted_rt'].tolist()])
		X = np.transpose(X)
		print(X.shape)

		reg = LinearRegression().fit(X, y)
		beta = reg.coef_[0]
		alpha = reg.intercept_
		print(beta, alpha)

		self.linear_beta = beta
		self.linear_alpha = alpha

		###Compute residuals###
		
		bench_PSMs_df['RT_residuals'] = bench_PSMs_df['predicted_rt'] - bench_PSMs_df['measured_rt']*beta - alpha
		ext_PSMs_df['RT_residuals'] = ext_PSMs_df['predicted_rt'] - ext_PSMs_df['measured_rt']*beta - alpha
		combined_PSMs_df['RT_residuals'] = combined_PSMs_df['predicted_rt'] - combined_PSMs_df['measured_rt']*beta - alpha
		
		self.confident_target_residuals = combined_PSMs_df.loc[(combined_PSMs_df['Decoy'] == 1) & (combined_PSMs_df['perc_qvalue'] <= 0.01)]['RT_residuals'].tolist()
		
		print('Confident targets: ', len(self.confident_target_residuals))

		confident_targets_residuals_std = np.std(np.array([self.confident_target_residuals]))
		confident_targets_residuals_mean = np.mean(np.array([self.confident_target_residuals]))

		print("Std :", confident_targets_residuals_std)
		print("Mean :", confident_targets_residuals_mean)

		self.pos_percentile_95 = confident_targets_residuals_mean + 1.645*confident_targets_residuals_std
		self.neg_percentile_95 = confident_targets_residuals_mean - 1.645*confident_targets_residuals_std
		self.pos_percentile_99 = confident_targets_residuals_mean + 2.326*confident_targets_residuals_std
		self.neg_percentile_99 = confident_targets_residuals_mean - 2.326*confident_targets_residuals_std
		
		return bench_PSMs_df, ext_PSMs_df		
			
			
			
			
			
