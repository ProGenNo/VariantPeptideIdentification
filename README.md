# Variant Peptide Identification
Code related to the "Retention time and fragmentation predictors increase confidence in variant peptide identification" publication: ...

### Requirements and Usage
Required software to reproduce the figures in the paper are **Snakemake** and **conda**, other required python libraries are specified at *diagnostic_plots_env.yaml* and the corresponding environment is automatically created by Snakemake.

* Download supplementary data from: ...
- Clone this repository
+ Adjust the paths to the data in Snakefile if necessary
- Create diagnostic plots by executing `snakemake --cores <# cores> --use-conda`

### Supplementary data
The supplementary data include the PSM and peptide level exports obtained from the analysis presented in this paper.
A proteomic dataset (provided by Wang et al. [1]) with MS/MS spectra from 3 samples was searched against 4 human protein databases (see publication for details):
* UniProt canonical proteome (UniProt)
- UniProt including isoforms (UniProt_isoforms)
+ Ensembl including isoforms (Ensembl_reference)
- Ensembl including isoforms and products of common genetic variants (extendedDB)

The exports contain a line per PSM/peptide with all relevant identifiers, features used for statistical evaluation, and confidence metrics.  
More detailed description of the columns of the exports are available at ...

The 4 used protein sequence databases are also included in the suppementary data.

### References
[1] D. Wang, B. Eraslan, T. Wieland, B. Hallström, T. Hopf, D. P. Zolg, J. Zecha, A. Asplund, L.-h. Li, C. Meng, M. Frejno, T. Schmidt, K. Schnatbaum, M. Wilhelm, F. Ponten, M. Uhlen, J. Gagneur, H. Hahne, and B. Kuster, “A deep proteome and transcriptome abundance atlas of 29 healthy human tissues,” Molecular Systems Biology, vol. 15, no. 2, p. e8503, 2019.
