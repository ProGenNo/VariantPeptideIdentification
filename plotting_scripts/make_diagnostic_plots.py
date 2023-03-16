import sys
import pandas as pd
from RT import *
from complex_plots import *
from peptide_plots import *

plots_directory_path = sys.argv[1]
database_name = sys.argv[2]

all_bench_PSMs_df = pd.read_csv(sys.argv[3])
all_bench_PSMs_df = all_bench_PSMs_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})

all_extra_PSMs_df = pd.read_csv(sys.argv[4])
all_extra_PSMs_df = all_extra_PSMs_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})

all_bench_Pept_df = pd.read_csv(sys.argv[5])
all_bench_Pept_df = all_bench_Pept_df.rename(columns={"benchmark_aggreg_score": "perc_score", "benchmark_aggreg_q-value": "perc_qvalue", "benchmark_aggreg_posterior_error_prob": "perc_pep"})

all_extra_Pept_df = pd.read_csv(sys.argv[6])
all_extra_Pept_df = all_extra_Pept_df.rename(columns={"extra_features_score": "perc_score", "extra_features_aggreg_q-value": "perc_qvalue", "extra_features_aggreg_posterior_error_prob": "perc_pep"})


RT = RT()

# RT residuals

all_bench_PSMs_df, all_extra_PSMs_df = RT.find_linear_model(all_bench_PSMs_df, all_extra_PSMs_df)

#########################################
######## 2D search space plots ##########
#########################################

ComplexPlots = SearchSpace2D(RT)

ComplexPlots.venn_PSMs(all_bench_PSMs_df, all_extra_PSMs_df, database_name, plots_directory_path, 'variant')
ComplexPlots.venn_PSMs(all_bench_PSMs_df, all_extra_PSMs_df, database_name, plots_directory_path, 'canonical')

all_bench_PSMs_df['RT_residuals'].mask(all_bench_PSMs_df['RT_residuals'] < -7000, -7000, inplace=True)
all_bench_PSMs_df['RT_residuals'].mask(all_bench_PSMs_df['RT_residuals'] > 7000, 7000, inplace=True)

all_extra_PSMs_df['RT_residuals'].mask(all_extra_PSMs_df['RT_residuals'] < -7000, -7000, inplace=True)
all_extra_PSMs_df['RT_residuals'].mask(all_extra_PSMs_df['RT_residuals'] > 7000, 7000, inplace=True)

ComplexPlots.all_TD_heatmaps(all_bench_PSMs_df, database_name, plots_directory_path, 'standard', 'spectra_angular_similarity', 'Spectra angular similarity', 300)

ComplexPlots.all_TD_psmType_heatmaps(all_bench_PSMs_df, database_name, plots_directory_path, 'standard', 'spectra_angular_similarity', 'Spectra angular similarity', 300, 'canonical')
ComplexPlots.all_TD_psmType_heatmaps(all_bench_PSMs_df, database_name, plots_directory_path, 'standard', 'spectra_angular_similarity', 'Spectra angular similarity', 300, 'variant')

ComplexPlots.conf_targets_acceptedGroups_RTresiduals_fragAng(all_bench_PSMs_df, all_extra_PSMs_df, database_name, plots_directory_path, ['canonical', 'variant'], 'standard', 'spectra_angular_similarity', 'Spectra angular similarity', 300)

ComplexPlots.conf_variants_scatter_RTresiduals_fragAng(all_bench_PSMs_df, all_extra_PSMs_df, database_name, plots_directory_path, ['variant'])

#########################################
############ Peptide plots ##############
#########################################

PeptidePlots = PeptidePlots()

PeptidePlots.venn_peptides(all_bench_Pept_df, all_extra_Pept_df, database_name, plots_directory_path, 'canonical')
PeptidePlots.venn_peptides(all_bench_Pept_df, all_extra_Pept_df, database_name, plots_directory_path, 'variant')

