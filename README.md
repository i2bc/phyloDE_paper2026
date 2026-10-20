Code to reproduce simulations and data analysis in 
"Interspecies Differential Gene Expression Analysis with Regularized Phylogenetic Linear Models" 
(Gallopin et al. 2026)

## Requirements

Packages:
```r
install.packages(c("here", "knitr", "ggplot2", "cowplot", "dplyr", "doParallel", "ape", "readxl", "data.table", "kableExtra", "readr"))
remotes::install_gitlab("sandve-lab/evemodel")
```

Last versions of `compcodeR` and `phyloDE`:
```r
remotes::install_github("pbastide/compcodeR", ref = "1f02368")
remotes::install_github("pbastide/phyloDE", ref = "40fa7f8")
```

## Content

* `mole_rat_RNAseq`: data and scripts

  * `00_data_format.R`: read and format the data
  
  * `01_plot_tree.R`: plot the tree transformed for several values of rho and lambda (Fig. 3).
  
  * `data`: raw and formatted data
    * `Daunesse2025`: raw data from [Daunesse et al. 2025](https://gitlab.com/maelle.daunesse/daunesse_et_al_2025)
    * `data_formatted`: formatted data, produced by `00_data_format.R`
  
  * `scripts`: simulation scripts
    * `01_simulations_compcodeR`: simulations using `compcodeR`
      * `01_simulations.R`: simulate interspecies count data using `compcodeR`
      * `02_analysis.R`: analyse the simulated data
      * `03_format_results.R`: read and format results
      * `04_plots.R`: plot figures 2, S1, S2, S3, S4, S5
        * **NOTE**: inputs needed to run the plot script are provided in the repository
        
    * `02_simulations_gaussian`: Gaussian simulation scripts
      * `01_data_sim.R`: simulate Gaussian data
      * `02_analysis.R`: analyse data
      * `03_plot.R`: plot figure 4
      * `04_test_one_simulation.R`: produce figures 1, S6, S7, S8
        * **NOTE**: inputs needed to run the plot script are provided in the repository
      
    * `03_simulations_trim`: Gaussian simulations to calibrate trim values
      * `01_gaussian_trim_test.R`: produce figure S9
      
    * `04_data_analysis`: re-analysis of the mole rat dataset
      * `01_data_analysis.R`: analyse the data with `phyloDE` and `evemodel`
      * `02_lfc_cre_vs_lfc_rnaseq.R`: CRE enrichment analysis
      * `03_plot_results_realdata_nbDE_nbDEconcor.R`: plot figure 5 and tables S1-6
    
  * `results`: scripts outputs
    * `2026-01-19_2026-01-19_simulations_results_mole_rat_RNAseq_data_Liver`: output from `mole_rat_RNAseq/scripts/01_simulations_compcodeR/01_simulations.R`
      * `simulation_parameters.RData`: saved simulation parameters
    * `2026-01-19_simulations_mole_rat_RNAseq_data_Liver`: output from `mole_rat_RNAseq/scripts/01_simulations_compcodeR/03_format_results.R`
      * `full_result_table_1_20.rds`: table of formatted results
    
    * `2026-01-15_gaussian_simulations`: output from `mole_rat_RNAseq/scripts/02_simulations_gaussian/01_data_sim.R`
      * `simulation_parameters.RData`: saved simulation parameters
    * `2026-01-15_2026-01-15_gaussian_simulations_results`: output from `mole_rat_RNAseq/scripts/02_simulations_gaussian/03_format_results.R`
      * `result_table.tsv`: table of formatted results
    
    * `2026-02-17_trim_simulations`: output from `mole_rat_RNAseq/scripts/03_simulations_trim/01_gaussian_trim_test.R`
      * `result_table.tsv`: table of formatted results
    
    * `2026-03-11_data_analysis` output from `mole_rat_RNAseq/scripts/04_data_analysis`
      * `mole_rat_RNAseq_data.tsv`: table of formatted results
      * `mole_rat_RNAseq_data_batch.tsv`: table of formatted results (with batch)
      * `fischer_lfc_cre_vs_lfc_rnaseq_by_tissue_branch.tsv`: table of CRE enrichment results
      * `mole_rat_RNAseq_TISSUE_DESIGN.tsv`: table of DE genes found by `phyloDE` for the TISSUE and DESIGN
      
* `crayfish_RNAseq`: data and scripts

  * `00_data_format.R`: read and format the data
  
  * `scripts`: simulation scripts
    * `01_simulations.R`: simulate interspecies count data using `compcodeR`
    * `02_analysis.R`: analyse the simulated data
    * `03_plots.R`: plot figure S10
      * **NOTE**: inputs needed to run the plot script are provided in the repository
    
  * `results`: scripts outputs
    * `stern2018_cpd`: formatted data, produced by `00_data_format.R`
    
    * `2026-01-19_simulations_crayfish`: output from `crayfish_RNAseq/scripts/01_simulations.R`
      * `params.RData`: saved simulation parameters
      * `designs.RData`: saved simulation parameters
      
    * `2026-01-19_2026-01-19_simulations_results_crayfish`: output from `crayfish_RNAseq/scripts/02_analysis.R`
      * `res.rds`: table of formatted results
      
* `mammals_RNAseq`: data and scripts

  * `001_data_format.R`: read and format the raw data from the remote repository
  * `002_data_format.R`: read and format the data
  
  * `data`: raw data as formatted by script `001_data_format.R`.
    
  * `scripts`: simulation scripts
    * `01_simulations.R`: simulate interspecies count data using `compcodeR`
    * `02_analysis.R`: analyse the simulated data
    * `03_plots.R`: plot figure S11
      * **NOTE**: inputs needed to run the plot script are provided in the repository
    
  * `results`: scripts outputs
    * `chen2019_primates_cpd`: formatted data, produced by `00_data_format.R`
    
    * `2026-02-09_simulations_mammals`: output from `mammals_RNAseq/scripts/01_simulations.R`
      * `params.RData`: saved simulation parameters
      * `designs.RData`: saved simulation parameters
      
    * `2026-02-09_2026-02-09_simulations_results_mammals`: output from `mammals_RNAseq/scripts/02_analysis.R`
      * `res.rds`: table of formatted results
    

