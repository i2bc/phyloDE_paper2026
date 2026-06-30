# Copyright (C) {2026} {PB, MG, AL}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## DESCRIPTION #################################################################
# Program to run differential expression analysis on simulated data
# NOTE: This is a computation intensive script
################################################################################

## INPUT #######################################################################
# Simulated data
# Saved in `mole_rat_RNAseq/results/YYYY-MM-DD_simulations_mole_rat_RNAseq_data_Liver`
################################################################################

## OUTPUT ######################################################################
# Analysis results, as `compcodeR` objects
# Saved in `mole_rat_RNAseq/results/YYYY-MM-DD_YYYY-MM-DD_simulations_results_mole_rat_RNAseq_data_Liver`
################################################################################

################################################################################
################################################################################

library(here)
library(compcodeR)
library(doParallel)

################################################################################
## Load Simulation Parameters
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"
data_name <- paste0("mole_rat_RNAseq_data_", tissue)
data_molerats <- readRDS(here(data_file, "data_formatted", data_name))

datestamp_day_simus <- "2026-01-19" ## Change here to match simulation date
simus_directory <- paste0(datestamp_day_simus, "_simulations_", data_name)
simus_directory <- file.path(here_dir, "mole_rat_RNAseq/results", simus_directory)
load(file.path(simus_directory, "simulation_parameters.RData"))

################################################################################
## File management
################################################################################
datestamp_day_anaysis <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day_anaysis, "_", datestamp_day_simus, "_simulations_results_", data_name)
results_directory <- file.path(here_dir, "mole_rat_RNAseq/results", results_directory)
dir.create(results_directory)


################################################################################
## Parallel Computations settings
################################################################################
## Iterations to analyze
Nmin <- 1
Nmax <- 10

## Required packages to pass to all nodes
reqpckg <- c("compcodeR", "here", "foreach", "phyloDE", "ape", "evemodel")

## Number of cores
Ncores <- 16

## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

################################################################################
## Loop on all settings
################################################################################

foreach (tree_type = all_tree_types) %:%
  foreach (cond_type = all_cond_types) %:%
  foreach (use_lengths = all_use_lengths) %:%
  foreach (model_process = c("BM", "OU")) %:%
  foreach (effect_size = all_effect_size) %:%
  foreach (prop_var_tree = all_prop_var_tree) %:%
  foreach (sel_strength = all_selection_strength) %:%
  foreach (fact_disp = all_fact_disp) %:%
  foreach (i = Nmin:Nmax) %:%
  foreach (method_type = c("limma", "phylolm", "phyloDE", "eve"), .packages = reqpckg, .verbose = TRUE) %dopar% {
    
    
    ## Remove some combinations of parameters
    if (prop_var_tree != sel_strength) return(NULL)
    
    ## Dataset names
    dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, sel_strength, fact_disp, sep = "_")
    
    dataset_info <- list(tree_type = tree_type,
                         cond_type = cond_type,
                         model_process = model_process,
                         use_lengths = use_lengths,
                         effect_size = effect_size,
                         prop_var_tree = prop_var_tree,
                         sel_strength = sel_strength,
                         fact_disp = fact_disp)
    
    ## Data file name
    data_file <- file.path(simus_directory, dataset)
    
    
    ####################################################################
    ## Limma cor
    ####################################################################
    if (method_type == "limma") {
      method <- "lengthNorm.limma"
      
      ## Length normalization when there are lengths only
      all_length_norm <- c("none")
      if (use_lengths == "with_lengths") all_length_norm <- c("TPM")
      ## Replicates correlation when there are replicates only
      all_blocks <- c("with_blocks")
      
      foreach (lnorm =  all_length_norm) %:%
        foreach (ltrans = c("log2")) %:%
        foreach (trend = c("no_trend", "with_trend")) %:%
        # foreach (trend = c("with_trend")) %:%
        foreach (block = all_blocks, .packages = reqpckg, .verbose = TRUE) %do% {
          
          ## Method id
          method_name <- paste(method, lnorm, ltrans, trend, block, sep = ".")
          
          method_info <- list(method = method,
                              lengthNormalization = lnorm,
                              transformation = ltrans,
                              methodProcess = NA,
                              eBayes = TRUE,
                              trend = trend,
                              cor = ifelse(block == "no_blocks", "no_cor", "cor"),
                              regCor = TRUE,
                              nulldisp = NA,
                              ddf_method = NA,
                              rep = i)
          method_info$methodMod <- paste0(method_info$method, "_",
                                          method_info$trend, "_",
                                          method_info$cor)
          
          ## If replicate correlation, find the right blocks
          is_block <- NULL
          if (block == "with_blocks") is_block <- "id.species"
          
          ## Sanity check : if result file already exists, do not run the analysis
          res_file_name <- file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds"))
          if (file.exists(res_file_name)) return(NULL)
          
          ## Run analysis
          el_time <- system.time(
            runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                       result.extent = method_name,
                       Rmdfunction = paste0(method, ".createRmd"),
                       output.directory = results_directory,
                       norm.method = "TMM",
                       length.normalization = lnorm,
                       data.transformation = ltrans,
                       trend = (trend == "with_trend"),
                       block.factor = is_block)
          )
          method_info$time <- el_time[3]
          
          ## Save parameters
          res_tmp <- readRDS(res_file_name)
          res_tmp@method.names <- c(res_tmp@method.names, dataset_info, method_info)
          saveRDS(res_tmp, file = res_file_name)
          
          ## For the first iteration only, generate the analysis code
          if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
        }
    }
    
    ##################################################################
    ## phylolm - with lengths
    ##################################################################
    if (method_type == "phylolm") {
      method <- "phylolm"
      
      if (tree_type != "no_tree") { ## Run phylolm only when there is a tree
        
        ## Length normalization when there are lengths only
        all_length_norm <- c("none")
        if (use_lengths == "with_lengths") all_length_norm <- c("TPM")
        
        foreach (sproc = c("BM", "OUfixedRoot")) %:%
          #foreach (sproc = c("OUfixedRoot")) %:%
          foreach (lnorm = all_length_norm) %:%
          foreach (ltrans = c("log2")) %:%
          foreach(i = Nmin:Nmax, .packages = reqpckg, .verbose = TRUE) %do% {
            
            ## Method id
            method_name <- paste(method, lnorm, ltrans, sproc, sep = ".")
            
            method_info <- list(method = method,
                                lengthNormalization = lnorm,
                                transformation = ltrans,
                                methodProcess = sproc,
                                eBayes = FALSE,
                                trend = NA,
                                cor = NA,
                                regCor = NA,
                                nulldisp = NA,
                                ddf_method = NA,
                                rep = i)
            method_info$methodMod <- paste0(method_info$method, "_",
                                            method_info$methodProcess)
            ## Sanity check : if result file already exists, do not run the analysis
            res_file_name <- file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds"))
            if (file.exists(res_file_name)) return(NULL)
            
            ## Run analysis
            el_time <- system.time(
              runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                         result.extent = method_name,
                         Rmdfunction = paste0(method, ".createRmd"),
                         output.directory = results_directory,
                         norm.method = "TMM",
                         model = sproc,
                         measurement_error = TRUE,
                         extra.design.covariates = NULL,
                         length.normalization = lnorm,
                         data.transformation = ltrans,
                         lower.bound = list(sigma2_error = (.Machine$double.eps)^0.5),
                         upper.bound = list(lambda = 1 / (1 + (.Machine$double.eps)^0.5)),
                         REML = TRUE)
            )
            method_info$time <- el_time[3]
            
            ## Save parameters
            res_tmp <- readRDS(res_file_name)
            res_tmp@method.names <- c(res_tmp@method.names, dataset_info, method_info)
            saveRDS(res_tmp, file = res_file_name)
            
            ## For the first iteration only, generate the analysis code
            if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
          }
      }
    }
    
    
    ##################################################################
    ## phyloDE - with lengths
    ##################################################################
    if (method_type == "phyloDE") {
      method <- "phyloDE"
      
      if (tree_type != "no_tree") { ## Run phyloDE only when there is a tree
        
        ## Length normalization when there are lengths only
        all_length_norm <- c("none")
        if (use_lengths == "with_lengths") all_length_norm <- c("TPM")
        
        foreach (sproc = c("BM", "OUfixedRoot")) %:%
          foreach (lnorm = all_length_norm) %:%
          foreach (ltrans = c("log2")) %:%
          foreach (regCor = c(TRUE, FALSE)) %:%
          foreach (eBayes = c(TRUE, FALSE)) %:%
          foreach (ddf_method = c("Samples")) %:%
          foreach (trend = c("no_trend", "with_trend"), .packages = reqpckg, .verbose = TRUE) %do% {
            
            if (trend == "no_trend" && eBayes == FALSE) return(NULL)
            
            ## Method id
            method_name <- paste(method, lnorm, ltrans, sproc, trend, paste0("eBayes_", eBayes), paste0("regCor_", regCor), paste0("ddf_method_", ddf_method), sep = ".")
            
            method_info <- list(method = method,
                                lengthNormalization = lnorm,
                                transformation = ltrans,
                                methodProcess = sproc,
                                eBayes = eBayes,
                                trend = trend,
                                cor = NA,
                                regCor = ifelse(regCor, "reg", "no_reg"),
                                nulldisp = NA,
                                ddf_method = ddf_method,
                                rep = i)
            method_info$methodMod <- paste0(method_info$method, "_",
                                            method_info$methodProcess, "_",
                                            paste0("eBayes_", method_info$eBayes), "_",
                                            paste0("ddf_method_", method_info$ddf_method), "_",
                                            method_info$trend, "_",
                                            method_info$regCor)
            
            ## Sanity check : if result file already exists, do not run the analysis
            res_file_name <- file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds"))
            if (file.exists(res_file_name)) return(NULL)
            
            ## Run analysis
            el_time <- system.time(
              runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                         result.extent = method_name,
                         Rmdfunction = paste0(method, ".createRmd"),
                         output.directory = results_directory,
                         norm.method = "TMM",
                         model = sproc,
                         measurement_error = TRUE,
                         extra.design.covariates = NULL,
                         length.normalization = lnorm,
                         data.transformation = ltrans,
                         use.eBayes = eBayes,
                         trend = (trend == "with_trend"),
                         regularize.correlation = regCor)
            )
            method_info$time <- el_time[3]
            
            ## Save parameters
            res_tmp <- readRDS(res_file_name)
            res_tmp@method.names <- c(res_tmp@method.names, dataset_info, method_info)
            saveRDS(res_tmp, file = res_file_name)
            
            ## For the first iteration only, generate the analysis code
            if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
          }
      }
    }
    
    ##################################################################
    ## evemodel- with lengths
    ##################################################################
    if (method_type == "eve") {
      method <- "evemodel"
      
      if (tree_type != "no_tree") { ## Run phyloDE only when there is a tree
        
        ## Length normalization when there are lengths only
        all_length_norm <- c("none")
        if (use_lengths == "with_lengths") all_length_norm <- c("TPM")
        
        foreach (lnorm = all_length_norm) %:%
          foreach (ltrans = c("log2")) %:%
          foreach (emp = c(TRUE), .packages = reqpckg, .verbose = TRUE) %do% {
            
            ## Method id
            method_name <- paste(method, ifelse(emp, "emp", "chi2"), lnorm, ltrans, sep = ".")
            
            method_info <- list(method = method,
                                lengthNormalization = lnorm,
                                transformation = ltrans,
                                methodProcess = "OU",
                                eBayes = FALSE,
                                trend = NA,
                                cor = NA,
                                regCor = NA,
                                nulldisp = ifelse(emp, "emp", "chi2"),
                                ddf_method = NA,
                                rep = i)
            method_info$methodMod <- paste0(method_info$method, "_",
                                            method_info$nulldisp)
            
            ## Sanity check : if result file already exists, do not run the analysis
            res_file_name <- file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds"))
            if (file.exists(res_file_name)) return(NULL)
            
            ## Run analysis
            el_time <- system.time(
              runDiffExp(data.file = paste0(data_file, "_", i, ".rds"),
                         result.extent = method_name,
                         Rmdfunction = "evemodel.twoThetaTest.createRmd",
                         output.directory = results_directory,
                         norm.method = "TMM",
                         length.normalization = "TPM",
                         empirical.p.values = emp,
                         n.genes.null.dist = 20000,
                         upperBound = c(theta = Inf, sigma2 = Inf, alpha = log(.Machine$double.xmax^0.8) / 2 / 1))
            )
            method_info$time <- el_time[3]
            
            ## Save parameters
            res_tmp <- readRDS(res_file_name)
            res_tmp@method.names <- c(res_tmp@method.names, dataset_info, method_info)
            saveRDS(res_tmp, file = res_file_name)
            
            ## For the first iteration only, generate the analysis code
            if (i == 1) generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", i, "_", method_name, ".rds")), results_directory)
          }
      }
    }
  }


## Stop the parallel cluster
stopCluster(cl)

