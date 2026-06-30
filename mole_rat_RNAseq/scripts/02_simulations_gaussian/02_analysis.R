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
# Program to analyse the Gaussian data
################################################################################

## INPUT #######################################################################
# Simulated data
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_gaussian_simulations`
################################################################################

## OUTPUT ######################################################################
# Results of the analyses
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_YYY-MM-DD_gaussian_simulations_results`
################################################################################

################################################################################
################################################################################

library(here)
library(compcodeR)
library(ape)
library(evemodel)
library(phyloDE)
library(doParallel)

################################################################################
## File management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"

datestamp_day_simu <- "2026-01-15"
simus_directory <- paste0(datestamp_day_simu, "_gaussian_simulations")
simus_directory <- file.path(here_dir, "mole_rat_RNAseq/results", simus_directory)
load(file.path(simus_directory, "simulation_parameters.RData"))

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day, "_", datestamp_day_simu, "_gaussian_simulations_results")
results_directory <- file.path(here_dir, "mole_rat_RNAseq/results", results_directory)
dir.create(results_directory, showWarnings = FALSE)

result_file <- file.path(results_directory, "result_table.tsv")

## result list to store methods already computed
result_list <- file.path(results_directory, "result_list.txt")
if (!file.exists(result_list)) file.create(result_list)

## tree without replicates
tree <- read.tree(file = file.path(here_dir, "mole_rat_RNAseq_data/species_tree_mya.nwk"))
tree$edge.length <- tree$edge.length / max(diag(vcv(tree)))

################################################################################
## Utility function
################################################################################
# Function to get parent branches of the set of tips.
# tips must be named 
get_parent_branches <- function(phy, tips) {
  # if only one tip, just get the parent branch
  if (length(tips) == 1) {
    return((phy$edge[, 2] == which(phy$tip.label == tips)) + 0)
  }
  # if several tips, get the MRCA, then descendants
  par <- ape::getMRCA(phy, tips)
  par <- c(par, phytools::getDescendants(phy, par))
  return(phy$edge[, 2] %in% par + 0)
}

getStats <- function(pvals, ntrue = ndiff, th = 0.05, name = "") {
  allTrue <- 1:length(pvals) <= ntrue
  TT <- sum(allTrue)
  
  TP <- sum(pvals <= th & allTrue)
  FP <- sum(pvals <= th & !allTrue)
  TN <- sum(pvals > th & !allTrue)
  FN <- sum(pvals > th & allTrue)
  
  N <- TP + FP + TN + FN
  S <- (TP + FN)/N
  P <- (TP + FP)/N
  
  if (FP + TP > 0) {
    FDR <- FP / (FP + TP)
    MCC <- (TP/N - S*P)/sqrt(P*S*(1 - S)*(1 - P))
  } else {
    FDR <- 0
    MCC <- 0
  }
  
  res <- list(
    FP = FP,
    TP = TP,
    TN = TN,
    FN = FN,
    FDR = FDR,
    TPR = TP / TT,
    MCC = MCC
  )
  names(res) <- paste0(names(res), name)
  
  return(res)
}

is_computed <- function(res_file_name) {
  if (isTRUE(any(grepl(res_file_name, readLines(result_list))))) {
    return(TRUE)
  }
  return(FALSE)
}

write_computed <- function(res_file_name) {
  result_list_con <- file(result_list, "a")
  writeLines(res_file_name, result_list_con)
  close(result_list_con)
}

write_table_line <- function(resultmethod, result_file) {
  if (file.exists(result_file)) {
    write.table(resultmethod, result_file, sep='\t', quote = FALSE, append = TRUE, col.names = FALSE)
  } else {
    write.table(resultmethod, result_file, sep='\t', quote = FALSE, append = TRUE, col.names = TRUE)
  }
}

################################################################################
## Parallel Computations settings
################################################################################
## Iterations to analyze
Nmin <- 1
Nmax <- 1

## Required packages to pass to all nodes
reqpckg <- c("compcodeR", "here", "foreach", "phyloDE", "evemodel")

## Number of cores
Ncores <- 5

## Register nodes
cl <- makeCluster(Ncores, outfile = "")
registerDoParallel(cl)

################################################################################
## Simulations BM
################################################################################

foreach (tree_type = all_trees) %:%
  foreach (branch = all_branch) %:%
  foreach (lambda_error_sim = all_lambda_error) %:%
  foreach (alpha_sim = all_alpha) %:%
  foreach (snr = all_snr) %:%
  foreach (irep = Nmin:Nmax) %:%
  foreach (method_name = c("limma", "phylolmBM", "phylolmOU", "phyloDEBM", "phyloDEOU", "eve")) %dopar% {
    
    if ((lambda_error_sim != lambda_error_base) && (alpha_sim != alpha_base)) return(NULL)
    
    ################################################################################
    ## select tree
    
    tree_sim <- switch(tree_type, original = tree_original, three_rep = tree_3rep)
    colDataTree <- colData[match(tree_sim$tip.label, rownames(colData)), ]
    
    ################################################################################
    ## dataset
    dataset <- paste(branch, tree_type, lambda_error_sim, alpha_sim, snr, irep, sep = "_")
    dataset_info <- data.frame(branch = branch,
                               tree_type = tree_type,
                               lambda_error_sim = lambda_error_sim,
                               alpha_sim = alpha_sim,
                               snr = snr,
                               irep = irep)
    
    ################################################################################
    ## load data
    sim_data <- read.table(file.path(simus_directory, paste0(dataset, ".tsv")))
    sim_data <- as.matrix(sim_data)
    
    ################################################################################
    ## design
    if (branch == "Hgla") {
      condition <- grepl("Hgla", colnames(sim_data)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Hgla"))
    } else if (branch == "MR") {
      condition <- grepl("Hgla|Fdam", colnames(sim_data)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Hgla", "Fdam"))
    }
    
    ## design linear model
    design_formula <- as.formula(~ condition)
    design_data <- data.frame(condition = condition)
    rownames(design_data) <- colnames(sim_data)
    design_data$condition <- factor(design_data$condition)
    design <- model.matrix(design_formula, design_data)
    
    ################################################################################
    ## limma
    if (method_name == "limma") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        block <- colDataTree$id.species
        corfit <- limma::duplicateCorrelation(sim_data, design = design, block = block, ndups = 1)
        rho <- corfit$consensus
        fitmethod <- limma::lmFit(sim_data, design = design,
                                  correlation = rho, block = block)
        fitmethod <- limma::eBayes(fitmethod, trend = FALSE)
        methodPvalues <- fitmethod$p.value[, ncol(fitmethod$p.value)]
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
    
    ################################################################################
    ## phyloDE - BM
    if (method_name == "phyloDEBM") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        fitmethod <- phyloDE::phylolmFit(sim_data,
                                         design = design,
                                         phy = tree,
                                         model = "BM",
                                         measurement_error = TRUE,
                                         use_consensus = TRUE,
                                         ddf_method = "Samples",
                                         trim = c(0.25, 0.05),
                                         ncores = 1)
        
        fitmethod <- limma::eBayes(fitmethod, trend = FALSE)
        methodPvalues <- fitmethod$p.value[, ncol(fitmethod$p.value)]
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
    
    ################################################################################
    ## phyloDE - OU
    if (method_name == "phyloDEOU") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        fitmethod <- phyloDE::phylolmFit(sim_data,
                                         design = design,
                                         phy = tree,
                                         model = "OUfixedRoot",
                                         measurement_error = TRUE,
                                         use_consensus = TRUE,
                                         ddf_method = "Samples",
                                         trim = c(0.25, 0.05),
                                         ncores = 1)
        
        fitmethod <- limma::eBayes(fitmethod, trend = FALSE)
        methodPvalues <- fitmethod$p.value[, ncol(fitmethod$p.value)]
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
    
    ################################################################################
    ## phylolm - BM
    if (method_name == "phylolmBM") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        fitmethod <- phyloDE::phylolmFit(sim_data,
                                         design = design,
                                         phy = tree,
                                         model = "BM",
                                         measurement_error = TRUE,
                                         use_consensus = FALSE,
                                         ddf_method = "Samples",
                                         ncores = 1)
        
        methodTvalues <- fitmethod$coef[, ncol(fitmethod$coef)] / fitmethod$stdev.unscaled[, ncol(fitmethod$stdev.unscaled)] / fitmethod$sigma
        methodPvalues <- 2 * pt(-abs(methodTvalues), df = fitmethod$df.residual)
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
    
    ################################################################################
    ## phylolm - OU
    if (method_name == "phylolmOU") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        fitmethod <- phyloDE::phylolmFit(sim_data,
                                         design = design,
                                         phy = tree,
                                         model = "OUfixedRoot",
                                         measurement_error = TRUE,
                                         use_consensus = FALSE,
                                         ddf_method = "Samples",
                                         ncores = 1)
        
        methodTvalues <- fitmethod$coef[, ncol(fitmethod$coef)] / fitmethod$stdev.unscaled[, ncol(fitmethod$stdev.unscaled)] / fitmethod$sigma
        methodPvalues <- 2 * pt(-abs(methodTvalues), df = fitmethod$df.residual)
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
    
    ################################################################################
    ## EVE
    if (method_name == "eve") {
      
      res_file_name <- file.path(results_directory, paste0(dataset, "_", irep, "_", method_name, ".tsv"))
      if (is_computed(res_file_name)) return(NULL)
      
      el_time <- system.time({
        ## first fit
        colSpecies <- sub("_.*$","", colnames(sim_data))
        fitEVE <- evemodel::twoThetaTest(tree = tree, gene.data = sim_data,
                                         isTheta2edge = isTheta2edge, colSpecies = colSpecies,
                                         cores = 1)
        ## median parameters of the null model estimated from the data
        params <- apply(fitEVE$oneThetaRes$par, 2, function(x) median(x, na.rm = TRUE))
        ## run null simulations
        set.seed(1289)
        data_sim_eve <- evemodel::simOneTheta(n = 20000, tree = tree, colSpecies = colSpecies,
                                              theta = params["theta"], sigma2 = params["sigma2"], 
                                              alpha = params["alpha"], beta = params["beta"])
        
        # run the evemodel to test for regulatory activity shift on the selected branch(es)
        fitEVEsim <- evemodel::twoThetaTest(tree = tree, gene.data = data_sim_eve,
                                            isTheta2edge = isTheta2edge, colSpecies = colSpecies,
                                            cores = 1)
        # Compute empirical p-values from simulation and write the TwoThetaTest results to file
        shift <- fitEVE$twoThetaRes$par[,"theta2"] - fitEVE$twoThetaRes$par[,"theta1"]
        methodPvalues <- qvalue::empPvals(fitEVE$LRT, fitEVEsim$LRT)
        methodPvaluesadj <- p.adjust(methodPvalues, method = 'BH')
      })
      
      resultmethod <- cbind(dataset_info, getStats(methodPvaluesadj), getStats(methodPvalues, name = ".raw"))
      resultmethod$time <- el_time[3]
      resultmethod$method <- method_name
      
      write_table_line(resultmethod, result_file)
      write_computed(res_file_name)
    }
  }

## Stop the parallel cluster
stopCluster(cl)

