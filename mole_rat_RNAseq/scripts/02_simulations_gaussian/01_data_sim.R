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
# Program to simulate gaussain data using the mole rats tree
################################################################################

## INPUT #######################################################################
# Formatted tree
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Simulated data
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_gaussian_simulations`
################################################################################

################################################################################
################################################################################

library(here)
library(compcodeR)
library(ape)

################################################################################
## Load Data
################################################################################
## read data
data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"
data_name <- paste0("mole_rat_RNAseq_data_", tissue)
data_molerats <- readRDS(here(data_file, "data_formatted", data_name))

## load tree and normalize
colData <- data_molerats@sample.annotations
tree_original <- data_molerats@tree
tree_original$edge.length <- tree_original$edge.length / max(diag(vcv(tree_original)))
t_original_tree <- phyloDE:::tree_height(tree_original)
# plot(tree_original)

## 3 replicates only
tree_3rep <- tree_original
tree_3rep <- drop.tip(tree_3rep, tree_original$tip.label[c(1:3, 10:11, 12:15, 19, 20:21, 25)])
# plot(tree_3rep)

all_trees <- c("original", "three_rep")

################################################################################
## Result files
################################################################################
here_dir <- here()

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day, "_gaussian_simulations")
results_directory <- file.path(here_dir, "mole_rat_RNAseq/results", results_directory)
dir.create(results_directory, showWarnings = FALSE)

################################################################################
## Utility function
################################################################################
get_lambda_error_OU <- function(phy, sigma2, sigma2_error, alpha) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                               parameters = list(alpha = alpha))$tree
  tilde_t <- phyloDE:::tree_height(tree_model) / (2 * alpha)
  lambda_ou_error <- phyloDE:::get_lambda_error(sigma2, sigma2_error, tilde_t)
  return(lambda_ou_error)
}

get_total_variance_BM <- function(sigma2, sigma2_error, t_tree) {
  return(sigma2 * t_tree + sigma2_error)
}

get_total_variance_OU <- function(phy, sigma2, sigma2_error, alpha) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                               parameters = list(alpha = alpha))$tree
  tilde_t <- phyloDE:::tree_height(tree_model) / (2 * alpha)
  return(get_total_variance_BM(sigma2, sigma2_error, tilde_t))
}

get_phylo_variance_OU <- function(phy, sigma2, alpha) {
  tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                               parameters = list(alpha = alpha))$tree
  tilde_t <- phyloDE:::tree_height(tree_model) / (2 * alpha)
  return(sigma2 * tilde_t)
}

get_parameters_BM <- function(lambda, t_tree, sigma2_total = 1) {
  return(list(sigma2_tree = lambda * sigma2_total / t_tree,
              sigma2_error = (1 - lambda) * sigma2_total))
}

get_parameters_OU <- function(alpha, lambda, phy, sigma2_total = 1) {
  if (alpha == 0) return(get_parameters_BM(lambda, phyloDE:::tree_height(phy), sigma2_total))
  var_ou <- get_phylo_variance_OU(phy, 1.0, alpha)
  return(list(sigma2_tree = lambda * sigma2_total / var_ou,
              sigma2_error = (1 - lambda) * sigma2_total))
}

################################################################################
## Parameters
################################################################################
nrep <- 50
ngenes <- 10000
all_branch <- c("MR", "Hgla")

all_lambda_error <- c(0.1, 0.25, 0.5, 0.75, 0.9)
lambda_error_base <- 0.75
all_alpha <- log(2) / c(0.05, 0.1, 0.5, 1, 5)
alpha_base <- log(2) / 0.5

## params OU
sapply(all_alpha, function(x) sapply(all_lambda_error, function(y) get_parameters_OU(x, y, tree_original)))

## differentially expressed genes
ndiff <- 150
all_snr <- c(2, 3, 4)

################################################################################
## Simulations
################################################################################

for (tree_type in all_trees) {
  tree_sim <- switch(tree_type, original = tree_original, three_rep = tree_3rep)
  colDataTree <- colData[match(tree_sim$tip.label, rownames(colData)), ]
  for (branch in all_branch) {
    ## Design
    if (branch == "Hgla") {
      condition <- colDataTree$conditionHgla
    } else if (branch == "MR") {
      condition <- colDataTree$conditionMR
    }
    for (lambda_error_sim in all_lambda_error) {
      for (alpha_sim in all_alpha) {
        if ((lambda_error_sim != lambda_error_base) && (alpha_sim != alpha_base)) next
        for (snr in all_snr) {
          params_sim <- get_parameters_OU(alpha_sim, lambda_error_sim, tree_sim)
          set.seed(1289)
          for (irep in 1:nrep) {
            
            ################################################################################
            ## Simulate OU
            sim_data <- t(phylolm::rTrait(ngenes, phy = tree_sim,
                                          model = "OU", parameters = list(sigma2 = params_sim$sigma2_tree, alpha = alpha_sim)))
            rownames(sim_data) <- paste0("g", 1:nrow(sim_data))
            
            ################################################################################
            ## Add noise
            sim_data <- sim_data + matrix(rnorm(prod(dim(sim_data)), 0, sqrt(params_sim$sigma2_error)), ncol = ncol(sim_data))
            
            ################################################################################
            ## Add effect size
            total_variance <- get_total_variance_OU(phy = tree_sim, sigma2 = params_sim$sigma2_tree, sigma2_error = params_sim$sigma2_error, alpha = alpha_sim)
            effect_size <- sqrt(total_variance) * snr
            sim_data[1:ndiff, ] <- sim_data[1:ndiff, ] + effect_size * matrix(rep(1, ndiff), ncol = 1) %*% matrix(as.numeric(condition) - 1, nrow = 1)
            
            ################################################################################
            ## Save
            dataset <- paste(branch, tree_type, lambda_error_sim, alpha_sim, snr, irep, sep = "_")
            write.table(sim_data, file = file.path(results_directory, paste0(dataset, ".tsv")))
          }
        }
      }
    }
  }
}

save(nrep,
     ngenes,
     all_branch,
     all_lambda_error,
     all_alpha,
     lambda_error_base,
     alpha_base,
     ndiff,
     all_snr,
     tree_original,
     tree_3rep,
     all_trees,
     colData,
     file = file.path(results_directory, "simulation_parameters.RData"))
