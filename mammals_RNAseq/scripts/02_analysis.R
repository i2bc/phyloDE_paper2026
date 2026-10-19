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
# Saved in `mammals_RNAseq/results/YYYY-MM-DD_simulations_mammals_RNAseq`
################################################################################

## OUTPUT ######################################################################
# Analysis results, as `compcodeR` objects
# Saved in `mammals_RNAseq/results/YYYY-MM-DD_YYYY-MM-DD_simulations_results_mammals`
################################################################################

library(here)
library(ape)
library(evemodel)
library(phylolm)
library(limma)
library(phyloDE)
library(parallel)

NB_THREADS <- 25

################################################################################
## Load Simulation Parameters
################################################################################
here_dir        <- here("mammals_RNAseq")
datestamp_simu  <- "2026-02-09" ## Change here to match simulation date
simus_directory <- file.path(here_dir, "results", paste0(datestamp_simu, "_simulations_mammals"))
params          <- readRDS(file.path(simus_directory, "params.rds"))
designs         <- readRDS(file.path(simus_directory, "designs.rds"))

################################################################################
## File management
################################################################################
datestamp_day     <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day, "_", datestamp_simu, "_simulations_results_mammals")
results_directory <- file.path(here_dir, "results", results_directory)
dir.create(results_directory)

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

getStats <- function(pvals, ntrue = ndiff, th = 0.05) {
  allTrue <- 1:length(pvals) <= ntrue
  TT <- sum(allTrue)
  
  TP <- sum(pvals <= th & allTrue)
  FP <- sum(pvals <= th & !allTrue)
  TN <- sum(pvals >  th & !allTrue)
  FN <- sum(pvals >  th & allTrue)
  
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
  
  res <- data.frame(
    FP  = FP,
    TP  = TP,
    TN  = TN,
    FN  = FN,
    FDR = FDR,
    TPR = TP / TT,
    MCC = MCC,
    th  = th
  )
  
  return(res)
}

remove_rep_from_tree <- function(tree) {
  species_names <- sub("\\..*", "", tree$tip.label)
  
  # Get one tip per species (first occurrence)
  unique_species <- unique(species_names)
  tips_to_keep <- sapply(unique_species, function(sp) {
    tree$tip.label[which(species_names == sp)[1]]
  })
  
  # Create species tree by keeping only one tip per species
  species_tree <- ape::keep.tip(tree, tips_to_keep)
  species_tree$tip.label <- sub("\\..*", "", species_tree$tip.label)
  
  species_tree
}

################################################################################
## Methods
################################################################################

methodFactory <- function(name) {
  switch(name, 
         phyloDEOU       = .phyloDEOU,
         phyloDEOUnotrim = .phyloDEOUnotrim,
         phyloDEBM       = .phyloDEBM,
         phyloDEBMnotrim = .phyloDEBMnotrim,
         eve             = .eve,
         limma           = .limma,
         phylolmBM       = .phylolmBM,
         phylolmOU       = .phylolmOU)
} 

.phyloDEOU <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "OUfixedRoot",
                                   measurement_error = TRUE,
                                   use_consensus     = TRUE,
                                   ddf_method        = "Samples",
                                   trim              = c(0.25, 0.05),
                                   ncores            = NB_THREADS)
  
  fitmethod <- limma::eBayes(fitmethod, trend = TRUE)
  pval      <- fitmethod$p.value[, ncol(fitmethod$p.value)]
  padj      <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.phyloDEOUnotrim <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "OUfixedRoot",
                                   measurement_error = TRUE,
                                   use_consensus     = TRUE,
                                   ddf_method        = "Samples",
                                   trim              = c(0, 1),
                                   ncores            = NB_THREADS)
  
  fitmethod <- limma::eBayes(fitmethod, trend = TRUE)
  pval      <- fitmethod$p.value[, ncol(fitmethod$p.value)]
  padj      <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.phyloDEBM <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "BM",
                                   measurement_error = TRUE,
                                   use_consensus     = TRUE,
                                   ddf_method        = "Samples",
                                   trim              = c(0.25, 0.05),
                                   ncores            = NB_THREADS)
  
  fitmethod <- limma::eBayes(fitmethod, trend = TRUE)
  pval      <- fitmethod$p.value[, ncol(fitmethod$p.value)]
  padj      <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.phyloDEBMnotrim <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "BM",
                                   measurement_error = TRUE,
                                   use_consensus     = TRUE,
                                   ddf_method        = "Samples",
                                   trim              = c(0, 1),
                                   ncores            = NB_THREADS)
  
  fitmethod <- limma::eBayes(fitmethod, trend = TRUE)
  pval      <- fitmethod$p.value[, ncol(fitmethod$p.value)]
  padj      <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.eve <- function(norm_cnts, design, tree) {
  
  norep_tree   <- remove_rep_from_tree(tree)
  colSpecies   <- sub("\\..*","", colnames(norm_cnts))
  isTheta2edge <- get_parent_branches(norep_tree, colSpecies[design[,2]==1])
  # first fit
  fitEVE <- evemodel::twoThetaTest(tree         = norep_tree, 
                                   gene.data    = norm_cnts,
                                   isTheta2edge = isTheta2edge,
                                   colSpecies   = colSpecies,
                                   cores        = NB_THREADS)
  # median parameters of the null model estimated from the data
  params <- apply(fitEVE$oneThetaRes$par, 2, function(x) median(x, na.rm = TRUE))
  
  # run null simulations
  set.seed(1289)
  data_sim_eve <- evemodel::simOneTheta(n          = 10000, 
                                        tree       = norep_tree, 
                                        colSpecies = colSpecies,
                                        theta      = params["theta"], 
                                        sigma2     = params["sigma2"], 
                                        alpha      = params["alpha"], 
                                        beta       = params["beta"])
  
  # run the evemodel to test for regulatory activity shift on the selected branch(es)
  fitEVEsim <- evemodel::twoThetaTest(tree         = norep_tree, 
                                      gene.data    = data_sim_eve,
                                      isTheta2edge = isTheta2edge, 
                                      colSpecies   = colSpecies,
                                      cores        = NB_THREADS)
  
  # Compute empirical p-values from simulation and write the TwoThetaTest results to file
  shift <- fitEVE$twoThetaRes$par[,"theta2"] - fitEVE$twoThetaRes$par[,"theta1"]
  pval  <- qvalue::empPvals(fitEVE$LRT, fitEVEsim$LRT)
  padj  <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.limma <- function(norm_cnts, design, tree) {
  block     <- sub("\\..*", "", rownames(design))
  corfit    <- limma::duplicateCorrelation(norm_cnts, design = design, block = block, ndups = 1)
  rho       <- corfit$consensus
  fitmethod <- limma::lmFit(norm_cnts, design = design, correlation = rho, block = block)
  fitmethod <- limma::eBayes(fitmethod, trend = TRUE)
  pval      <- fitmethod$p.value[, ncol(fitmethod$p.value)]
  padj      <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

.phylolmOU <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "OUfixedRoot",
                                   measurement_error = TRUE,
                                   use_consensus     = FALSE,
                                   ddf_method        = "Samples",
                                   ncores            = NB_THREADS)
  
  methodTvalues <- fitmethod$coef[, ncol(fitmethod$coef)] / fitmethod$stdev.unscaled[, ncol(fitmethod$stdev.unscaled)] / fitmethod$sigma
  pval          <- 2 * pt(-abs(methodTvalues), df = fitmethod$df.residual)
  padj          <- p.adjust(pval, method = 'BH')
  
  return(padj)
}
.phylolmBM <- function(norm_cnts, design, tree) {
  fitmethod <- phyloDE::phylolmFit(norm_cnts,
                                   design            = design,
                                   phy               = tree,
                                   model             = "BM",
                                   measurement_error = TRUE,
                                   use_consensus     = FALSE,
                                   ddf_method        = "Samples",
                                   ncores            = NB_THREADS)
  
  methodTvalues <- fitmethod$coef[, ncol(fitmethod$coef)] / fitmethod$stdev.unscaled[, ncol(fitmethod$stdev.unscaled)] / fitmethod$sigma
  pval          <- 2 * pt(-abs(methodTvalues), df = fitmethod$df.residual)
  padj          <- p.adjust(pval, method = 'BH')
  
  return(padj)
}

################################################################################
## Main Simulation Analysis
################################################################################
## For each combination of:
##   - Simulation parameters (effect size, model process, etc.)
##   - Species selection design (k species per condition)
##   - DE detection method (phyloDE, limma, eve, phylolm)
## Apply the method and compute performance statistics (FDR, TPR, MCC)
################################################################################

# Create expanded parameter grid combining all factors
expand_params <- expand.grid(
  id      = 1:nrow(params),
  method  = c("phyloDEOU", "phyloDEBM", "eve", "limma", "phylolmBM", "phylolmOU")
)

# Merge with simulation parameters
params$id     <- 1:nrow(params)
expand_params <- merge(expand_params, params)

res <- do.call(rbind, parallel::mclapply(1:nrow(expand_params), function(idx){
  
  # Extract parameters for this iteration
  method       <- as.character(expand_params$method[idx])
  dataset_file <- paste0(expand_params$dataset[idx], "_", expand_params$rep[idx], ".rds")
  
  # Load simulated dataset
  dataset   <- readRDS(file.path(simus_directory, dataset_file))
  tips_info <- dataset@sample.annotations
  norm_cnts <- log2(sweep(dataset@count.matrix, 2, tips_info$depth.factor, "/") + 0.5)
  
  # Create design matrix for linear model
  design_formula <- as.formula(~ condition)
  design_mat     <- model.matrix(design_formula, tips_info)
  if (any(design_mat[,2] == 2)) design_mat[,2] <- design_mat[,2] -1
  
  # Apply DE detection method and get adjusted p-values
  padj  <- methodFactory(method)(norm_cnts, design_mat, dataset@tree)
  ndiff <- sum(dataset@variable.annotations$differential.expression == 1)
  
  # Compute performance statistics
  result <- do.call(rbind, lapply(c(0.05, 0.1), function(th) {
    stats <- getStats(padj, ndiff, th = th)
    cbind(expand_params[idx, ], stats)
  }))
  
  return(result)
  
}, mc.cores = 1))

# Save results
saveRDS(res, file.path(results_directory, "res.rds"))
