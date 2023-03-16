DATABASES = ['extendedDB', 'UniProt', 'UniProt_isoforms', 'Ensembl_reference']
DATABASE = 'extendedDB' #or 'UniProt', 'UniProt_isoforms', 'Ensembl_reference'

rule all:
	input:
		in1=expand("Plots/{database}/Peptide_plots/confident_canonical_peptides_venn.png", database=DATABASE),
		in2="Plots/Experiments/experiments_numOfHits_Percolator_extendedDB.png",
		in3="Plots/Experiments/Peptide_plots/experiments_numOfHits_engine_score.png",
		in4="Plots/Combined_databases/RT_residuals_quants_jointplots.png",
		in5="Plots/Combined_databases/Peptide_plots/RT_residuals_quants_jointplots.png"
		

rule diagnostic_plots:
	input:
		in1="data/{database}/concatenated_ALL_benchmark_PSM_full_export",
		in2="data/{database}/concatenated_ALL_extra_PSM_full_export",
		in3="data/{database}/concatenated_ALL_benchmark_Peptide_full_export",
		in4="data/{database}/concatenated_ALL_extra_Peptide_full_export"
	output:
		out="Plots/{database}/Peptide_plots/confident_canonical_peptides_venn.png"
	params:
		plots_directory_path = "Plots/{database}",
		database_name = "{database}"
	conda:  "diagnostic_plots_env.yaml"
	shell:
                "mkdir -p {params.plots_directory_path}/2D_plots ; "
                "mkdir -p {params.plots_directory_path}/Peptide_plots ; "
                "python3 make_diagnostic_plots.py {params.plots_directory_path} {params.database_name} {input.in1} {input.in2} {input.in3} {input.in4}"
                
rule combined_databases_plots:
	input:
		in1="data/extendedDB/concatenated_ALL_benchmark_PSM_full_export",
		in2="data/extendedDB/concatenated_ALL_extra_PSM_full_export",
		in3="data/UniProt/concatenated_ALL_benchmark_PSM_full_export",
		in4="data/UniProt/concatenated_ALL_extra_PSM_full_export",
		in5="data/UniProt_isoforms/concatenated_ALL_benchmark_PSM_full_export",
		in6="data/UniProt_isoforms/concatenated_ALL_extra_PSM_full_export",
		in7="data/Ensembl_reference/concatenated_ALL_benchmark_PSM_full_export",
		in8="data/Ensembl_reference/concatenated_ALL_extra_PSM_full_export",
	output:
		out1="Plots/Combined_databases/RT_residuals_databases_TD_violin.png",
		out2="Plots/Combined_databases/spectra_angular_similarity_databases_TD_violin.png",
		out3="Plots/Combined_databases/RT_residuals_quants_jointplots.png"
	params:
		plots_directory_path = "Plots/Combined_databases"
	conda:  "diagnostic_plots_env.yaml"
	threads: 4
	shell:
		"mkdir -p ./{params.plots_directory_path} ; "
		"python3 Combine_databases.py {input.in1} {input.in2} {input.in3} {input.in4} {input.in5} {input.in6} {input.in7} {input.in8} {params.plots_directory_path}"

rule combined_databases_peptide_plots:
	input:
		in1="data/extendedDB/concatenated_ALL_benchmark_Peptide_full_export",
		in2="data/extendedDB/concatenated_ALL_extra_Peptide_full_export",
		in3="data/UniProt/concatenated_ALL_benchmark_Peptide_full_export",
		in4="data/UniProt/concatenated_ALL_extra_Peptide_full_export",
		in5="data/UniProt_isoforms/concatenated_ALL_benchmark_Peptide_full_export",
		in6="data/UniProt_isoforms/concatenated_ALL_extra_Peptide_full_export",
		in7="data/Ensembl_reference/concatenated_ALL_benchmark_Peptide_full_export",
		in8="data/Ensembl_reference/concatenated_ALL_extra_Peptide_full_export",
	output:
		out1="Plots/Combined_databases/Peptide_plots/RT_residuals_quants_jointplots.png"
	params:
		plots_directory_path = "Plots/Combined_databases/Peptide_plots"
	conda:  "diagnostic_plots_env.yaml"
	threads: 4
	shell:
		"mkdir -p ./{params.plots_directory_path} ; "
		"python3 Combine_databases.py {input.in1} {input.in2} {input.in3} {input.in4} {input.in5} {input.in6} {input.in7} {input.in8} {params.plots_directory_path}"
		
rule experiment_PSM_plots:
	input:
		in1="data/extendedDB/experiments/experiment_P010694_ALL_benchmark_PSM_export",
		in2="data/extendedDB/experiments/experiment_P010694_ALL_extra_PSM_export",
		in3="data/extendedDB/experiments/experiment_P010747_ALL_benchmark_PSM_export",
		in4="data/extendedDB/experiments/experiment_P010747_ALL_extra_PSM_export",
		in5="data/extendedDB/experiments/experiment_P013107_ALL_benchmark_PSM_export",
		in6="data/extendedDB/experiments/experiment_P013107_ALL_extra_PSM_export",
		in7="data/UniProt/experiments/experiment_P010694_ALL_benchmark_PSM_export",
		in8="data/UniProt/experiments/experiment_P010694_ALL_extra_PSM_export",
		in9="data/UniProt/experiments/experiment_P010747_ALL_benchmark_PSM_export",
		in10="data/UniProt/experiments/experiment_P010747_ALL_extra_PSM_export",
		in11="data/UniProt/experiments/experiment_P013107_ALL_benchmark_PSM_export",
		in12="data/UniProt/experiments/experiment_P013107_ALL_extra_PSM_export",
		in13="data/Ensembl_reference/experiments/experiment_P010694_ALL_benchmark_PSM_export",
		in14="data/Ensembl_reference/experiments/experiment_P010694_ALL_extra_PSM_export",
		in15="data/Ensembl_reference/experiments/experiment_P010747_ALL_benchmark_PSM_export",
		in16="data/Ensembl_reference/experiments/experiment_P010747_ALL_extra_PSM_export",
		in17="data/Ensembl_reference/experiments/experiment_P013107_ALL_benchmark_PSM_export",
		in18="data/Ensembl_reference/experiments/experiment_P013107_ALL_extra_PSM_export",
		in19="data/UniProt_isoforms/experiments/experiment_P010694_ALL_benchmark_PSM_export",
		in20="data/UniProt_isoforms/experiments/experiment_P010694_ALL_extra_PSM_export",
		in21="data/UniProt_isoforms/experiments/experiment_P010747_ALL_benchmark_PSM_export",
		in22="data/UniProt_isoforms/experiments/experiment_P010747_ALL_extra_PSM_export",
		in23="data/UniProt_isoforms/experiments/experiment_P013107_ALL_benchmark_PSM_export",
		in24="data/UniProt_isoforms/experiments/experiment_P013107_ALL_extra_PSM_export"
	output:
		out1="Plots/Experiments/experiments_numOfHits_Percolator_extendedDB.png",
		out2="Plots/Experiments/experiments_numOfHits_engine_score.png"
	params:
		plots_directory_path = "Plots/Experiments",
		analysis_level = "PSM"
	conda:  "diagnostic_plots_env.yaml"
	threads: 4
	shell:
		"mkdir -p ./{params.plots_directory_path} ; "
		"python3 experiment_plots.py {input.in1} {input.in2} {input.in3} {input.in4} {input.in5} {input.in6} {input.in7} {input.in8} {input.in9} {input.in10} {input.in11} {input.in12} {input.in13} {input.in14} {input.in15} {input.in16} {input.in17} {input.in18} {input.in19} {input.in20} {input.in21} {input.in22} {input.in23} {input.in24} {params.plots_directory_path} {params.analysis_level}"
		
rule experiment_peptide_plots:
	input:
		in1="data/extendedDB/experiments/experiment_P010694_ALL_benchmark_Peptide_export",
		in2="data/extendedDB/experiments/experiment_P010694_ALL_extra_Peptide_export",
		in3="data/extendedDB/experiments/experiment_P010747_ALL_benchmark_Peptide_export",
		in4="data/extendedDB/experiments/experiment_P010747_ALL_extra_Peptide_export",
		in5="data/extendedDB/experiments/experiment_P013107_ALL_benchmark_Peptide_export",
		in6="data/extendedDB/experiments/experiment_P013107_ALL_extra_Peptide_export",
		in7="data/UniProt/experiments/experiment_P010694_ALL_benchmark_Peptide_export",
		in8="data/UniProt/experiments/experiment_P010694_ALL_extra_Peptide_export",
		in9="data/UniProt/experiments/experiment_P010747_ALL_benchmark_Peptide_export",
		in10="data/UniProt/experiments/experiment_P010747_ALL_extra_Peptide_export",
		in11="data/UniProt/experiments/experiment_P013107_ALL_benchmark_Peptide_export",
		in12="data/UniProt/experiments/experiment_P013107_ALL_extra_Peptide_export",
		in13="data/Ensembl_reference/experiments/experiment_P010694_ALL_benchmark_Peptide_export",
		in14="data/Ensembl_reference/experiments/experiment_P010694_ALL_extra_Peptide_export",
		in15="data/Ensembl_reference/experiments/experiment_P010747_ALL_benchmark_Peptide_export",
		in16="data/Ensembl_reference/experiments/experiment_P010747_ALL_extra_Peptide_export",
		in17="data/Ensembl_reference/experiments/experiment_P013107_ALL_benchmark_Peptide_export",
		in18="data/Ensembl_reference/experiments/experiment_P013107_ALL_extra_Peptide_export",
		in19="data/UniProt_isoforms/experiments/experiment_P010694_ALL_benchmark_Peptide_export",
		in20="data/UniProt_isoforms/experiments/experiment_P010694_ALL_extra_Peptide_export",
		in21="data/UniProt_isoforms/experiments/experiment_P010747_ALL_benchmark_Peptide_export",
		in22="data/UniProt_isoforms/experiments/experiment_P010747_ALL_extra_Peptide_export",
		in23="data/UniProt_isoforms/experiments/experiment_P013107_ALL_benchmark_Peptide_export",
		in24="data/UniProt_isoforms/experiments/experiment_P013107_ALL_extra_Peptide_export"
	output:
		out1="Plots/Experiments/Peptide_plots/experiments_numOfHits_engine_score.png",
	params:
		plots_directory_path = "Plots/Experiments/Peptide_plots",
		analysis_level = "peptide"
	conda:  "diagnostic_plots_env.yaml"
	threads: 4
	shell:
		"mkdir -p ./{params.plots_directory_path} ; "
		"python3 experiment_plots.py {input.in1} {input.in2} {input.in3} {input.in4} {input.in5} {input.in6} {input.in7} {input.in8} {input.in9} {input.in10} {input.in11} {input.in12} {input.in13} {input.in14} {input.in15} {input.in16} {input.in17} {input.in18} {input.in19} {input.in20} {input.in21} {input.in22} {input.in23} {input.in24} {params.plots_directory_path} {params.analysis_level}"
