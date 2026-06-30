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
# Program to simulate data from empirical parameters
# using the mole rats dataset
################################################################################

## INPUT #######################################################################
# Formatted data
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Simulated data
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_simulations_mole_rat_RNAseq_data_Liver`
################################################################################

################################################################################
################################################################################

library(here)
library(compcodeR)
library(knitr)
library(ggplot2)
library(ape)

################################################################################
## Load Data
################################################################################
## read data
data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"
data_name <- paste0("mole_rat_RNAseq_data_", tissue)
data_molerats <- readRDS(here(data_file, "data_formatted", data_name))

counts_noNA <- data_molerats@count.matrix
length_noNA <- data_molerats@length.matrix
colData <- data_molerats@sample.annotations
tree_original <- data_molerats@tree

################################################################################
## Files management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
simus_directory <- paste0(datestamp_day, "_simulations_", data_name)
simus_directory <- file.path(here_dir, "mole_rat_RNAseq/results", simus_directory)
dir.create(simus_directory, showWarnings = FALSE)

################################################################################
## Analyse Data to get parameters
################################################################################
## remove zeros
small_counts <- rowSums(counts_noNA) <= 100
sum(small_counts) / nrow(counts_noNA)
counts_noNA <- counts_noNA[!small_counts, ]
length_noNA <- length_noNA[!small_counts, ]

## run DESeq2 to get parameters estimation
library(DESeq2)
# Correct by families
dds <- DESeqDataSetFromMatrix(counts_noNA, colData, ~ condition )
# normalisation factors for library size
dds <-  estimateSizeFactors(dds)
size_fac <- sizeFactors(dds)
# compute normalization with length
mat_size_fac <- matrix(size_fac, ncol=length(size_fac), nrow=length(counts_noNA[,1]), byrow=T)
normFactors <- (mat_size_fac*length_noNA) / exp(rowMeans(log(mat_size_fac*length_noNA)))
normalizationFactors(dds) <- as.matrix(normFactors)
# analysis
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
dds <- DESeq2::replaceOutliersWithTrimmedMean(dds)
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
res <- DESeq2::results(dds, independentFiltering = TRUE, cooksCutoff = TRUE)

################################################################################
## Get empirical moments for simulation
################################################################################
dispersions_molerats <- DESeq2::dispersions(dds)
seqdepth_molerats <- mean(colSums(counts_noNA))
minfact_molerats <- min(colSums(counts_noNA)/median(colSums(counts_noNA)))
maxfact_molerats <- max(colSums(counts_noNA)/median(colSums(counts_noNA)))
relmeans_molerats <- rowMeans(counts_noNA[,colData$condition == 1])
n.vars <- dim(counts_noNA)[1]

## plot mean var graph
plot(log(relmeans_molerats), log(dispersions_molerats))
points(log(relmeans_molerats)[1:150], log(dispersions_molerats)[1:150], col = "red")
hist(log(dispersions_molerats))
hist(log(dispersions_molerats)[1:150])

## replace NA in dispersion_molerats with values
anyNA(dispersions_molerats)

################################################################################
## Format Tree
################################################################################
n <- length(unique(colData$species)) # number of species
N <- length(tree_original$tip.label) # total number of observations (with replicates)

tree_original$edge.length <- tree_original$edge.length / max(diag(vcv(tree_original))[1:N]) # normalize tree to unit height

## 3 replicates only
tree_3rep <- tree_original
tree_3rep <- drop.tip(tree_3rep, tree_original$tip.label[c(1:3, 10:11, 12:15, 19, 20:21, 25)])
# plot(tree_3rep)

all_tree_types <- c("original", "three_rep")

################################################################################
## Format Lengths and compute empirical stats
################################################################################
species_names <- colData$species
unique_length <- length_noNA[, !duplicated(species_names)]

mean_lengths <- rowMeans(unique_length)
var_lengths <- apply(unique_length,1,var)
disp_lengths <- (var_lengths - mean_lengths) / mean_lengths^2 # size = mu^2/(var-mu), disp = 1 / size
disp_lengths[disp_lengths <= 0] <- runif(sum(disp_lengths <= 0), 0, quantile(disp_lengths, 0.01))

plot(log(mean_lengths), log(var_lengths))
abline(0,1)

################################################################################
## Empirical fraction of data explained by the tree
################################################################################
# Normalisation TPM
nf <- edgeR::calcNormFactors(counts_noNA / length_noNA, method = 'TMM')
lib.size <- colSums(counts_noNA / length_noNA) * nf
data.norm <- sweep((counts_noNA + 0.5) / length_noNA, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6
# Transformation log2
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(counts_noNA)
# phylolm pagel lambda analysis
all_fits_lambda <- apply(data.trans, 1, function(x) phylolm::phylolm(x ~ colData$condition, phy = tree_original, model = "lambda", upper.bound = list(lambda = 1 / (1 + (.Machine$double.eps)^0.5))))
lambdas <- sapply(all_fits_lambda, function(x) x$optpar)
all_means <- sapply(all_fits_lambda, function(x) x$coefficients[1])
hist(lambdas, breaks = 100)

lambda_median <- median(lambdas)
lambda_median

v_emp <- c(0.05^2, 0.1^2)
alpha_emp <- lambda_median * (lambda_median * (1 - lambda_median) / v_emp - 1)
beta_emp <- (1 - lambda_median) * (lambda_median * (1 - lambda_median) / v_emp - 1)
set.seed(1289)
lambdas_sim_low <- rbeta(length(lambdas), alpha_emp[1], beta_emp[1])
hist(lambdas_sim_low)
lambdas_sim_high <- rbeta(length(lambdas), alpha_emp[2], beta_emp[2])
hist(lambdas_sim_high)

# OU analysis to get empirical selection strength
oufits <- apply(data.trans, 1, function(x) phylolm::phylolm(x ~ colData$condition, phy = tree_original,
                                                            model = "OUfixedRoot", measurement_error = TRUE,
                                                            lower.bound = list(sigma2_error = (.Machine$double.eps)^0.5)))

get_lambda_error <- function(sigma2, sigma2_error, h_tree) {
  return(sigma2 * h_tree / (sigma2_error + sigma2 * h_tree))
}
get_gamma <- function(phyfit) {
  return(phyfit$sigma2 / 2 / phyfit$optpar)
}
tree_height <- function(tree) {
  return(max(ape::node.depth.edgelength(tree)))
}
get_lambda_error_OU <- function(alpha, sigma2, sigma2_error) {
  tree_model <- phylolm::transf.branch.lengths(tree_original, "OUfixedRoot",
                                               parameters = list(alpha = alpha))$tree
  tilde_t <- tree_height(tree_model) / (2 * alpha)
  lambda_ou_error <- get_lambda_error(sigma2, sigma2_error, tilde_t)
  return(lambda_ou_error)
}
get_lambda_error_OU_lm <- function(phyfit) {
  get_lambda_error_OU(phyfit$optpar, phyfit$sigma2, phyfit$sigma2_error)
}
alphas <- sapply(oufits, function(x) x$optpar)
# rhos <- get_rho(alphas, tree_height(tree))
sig2_err <- sapply(oufits, function(x) x$sigma2_error)
# gammas <- sapply(oufits, get_gamma)
lambdas_OU <- sapply(oufits, get_lambda_error_OU_lm)
hist(alphas, breaks = 100)
# hist(rhos, breaks = 100)
hist(log(2) / alphas, breaks = 100)
# hist(gammas, breaks = 100)
hist(lambdas_OU, breaks = 100)

alpha_median <- median(alphas)
alpha_median
log(2) / alpha_median
lambda_OU_median <- median(lambdas_OU)
lambda_OU_median

v_emp <- c(0.05^2, 0.1^2)
alpha_emp <- lambda_OU_median * (lambda_OU_median * (1 - lambda_OU_median) / v_emp - 1)
beta_emp <- (1 - lambda_OU_median) * (lambda_OU_median * (1 - lambda_OU_median) / v_emp - 1)
set.seed(1289)
lambdas_OU_sim_low <- rbeta(length(lambdas_OU), alpha_emp[1], beta_emp[1])
hist(lambdas_OU_sim_low)
lambdas_OU_sim_high <- rbeta(length(lambdas_OU), alpha_emp[2], beta_emp[2])
hist(lambdas_OU_sim_high)

v_emp <- c(1^2, 5^2)
theta_emp <- v_emp / alpha_median
alpha_emp <- alpha_median / theta_emp
set.seed(1289)
alphas_OU_sim_low <- rgamma(length(alphas), shape = alpha_emp[1], scale = theta_emp[1])
hist(alphas_OU_sim_low)
hist(log(2) / alphas_OU_sim_low)
alphas_OU_sim_high <- rgamma(length(alphas), shape = alpha_emp[2], scale = theta_emp[2])
hist(alphas_OU_sim_high)
hist(log(2) / alphas_OU_sim_high)

################################################################################
## Parameters for the simulation
################################################################################
n.diffexp <- 150

Nrep <- 50 # number of replicates
all_effect_size <- c(3.0, 4.0, 5.0)

all_selection_strength <- c("simlow", "simhigh", "median")

all_prop_var_tree <- c("simlow", "simhigh", "median")

all_fact_disp <- c(1)

all_use_lengths <- c("with_lengths")

all_cond_types <- c("MR", "Hgla")

################################################################################
## Simulations 
################################################################################

for (tree_type in all_tree_types) {
  for (cond_type in all_cond_types) {
    for (use_lengths in all_use_lengths) {
      for (model_process in c("BM", "OU")) {
        for (effect_size in all_effect_size) {
          for(prop_var_tree in all_prop_var_tree) {
            for(sel_strength in all_selection_strength) {
              for(fact_disp in all_fact_disp) {
                
                ## Choose tree
                tree_sim <- switch(tree_type, original = tree_original, three_rep = tree_3rep)
                colDataTree <- colData[match(tree_sim$tip.label, rownames(colData)), ]
                
                ## species id
                id_species <- colDataTree$species
                id_species <- factor(id_species)
                names(id_species) <- tree_sim$tip.label
                
                ## condition
                if (cond_type == "Hgla") {
                  id_cond <- colDataTree$conditionHgla
                } else if (cond_type == "MR") {
                  id_cond <- colDataTree$conditionMR
                }
                levels(id_cond) <- c(1, 2)
                names(id_cond) <- rownames(colDataTree)
                # reorder
                id_cond <- id_cond[match(tree_sim$tip.label, names(id_cond))]
                # plot(tree_sim)
                # tiplabels(pch = 21, col = id_cond, bg = id_cond)
                
                ## dataset
                dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, sel_strength, fact_disp, sep = "_")
                
                ## parameters
                pvt <- prop_var_tree
                sl <- sel_strength
                if (pvt != sl) break
                if (pvt == "simlow") {
                  if (model_process == "BM") pvt <- lambdas_sim_low
                  if (model_process == "OU") pvt <- lambdas_OU_sim_low
                } else if (pvt == "simhigh") {
                  if (model_process == "BM") pvt <- lambdas_sim_high
                  if (model_process == "OU") pvt <- lambdas_OU_sim_high
                } else if (pvt == "median") {
                  if (model_process == "BM") pvt <- lambda_median
                  if (model_process == "OU") pvt <- lambda_OU_median
                }
                if (model_process == "OU") {
                  if (sl == "simlow") {
                    sl <- alphas_OU_sim_low
                  } else if (sl == "simhigh") {
                    sl <- alphas_OU_sim_high
                  } else if (sl == "median") {
                    sl <- alpha_median
                  }
                }
                if (model_process == "BM") sl <- 0.0
                
                set.seed(18570823)
                for (i in 1:Nrep) {
                  simus_obj <- generateSyntheticData(dataset = dataset,
                                                     n.vars = n.vars,
                                                     n.diffexp = n.diffexp,
                                                     repl.id = i,
                                                     seqdepth = seqdepth_molerats,
                                                     minfact = minfact_molerats,
                                                     maxfact = maxfact_molerats,
                                                     relmeans = relmeans_molerats,
                                                     dispersions = dispersions_molerats * fact_disp,
                                                     fraction.upregulated = 0.5,
                                                     between.group.diffdisp = FALSE,
                                                     filter.threshold.total = 1,
                                                     filter.threshold.mediancpm = 0,
                                                     fraction.non.overdispersed = 0,
                                                     random.outlier.high.prob = 0,
                                                     random.outlier.low.prob = 0,
                                                     single.outlier.high.prob = 0,
                                                     single.outlier.low.prob = 0,
                                                     effect.size = effect_size,
                                                     output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                                     tree = tree_sim,
                                                     prop.var.tree = pvt,
                                                     id.condition = id_cond,
                                                     id.species = id_species,
                                                     model.process = model_process,
                                                     selection.strength = sl,
                                                     use.lengths = use_lengths == "with_lengths",
                                                     lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                                     lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL
                  )
                  
                  if (i == 1) { ## Show summary only for the first iteration
                    summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                              output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

save(Nrep,
     all_effect_size,
     all_use_lengths,
     all_tree_types,
     all_cond_types,
     all_prop_var_tree,
     all_selection_strength,
     all_fact_disp,
     dispersions_molerats,
     relmeans_molerats,
     seqdepth_molerats,
     minfact_molerats,
     maxfact_molerats,
     mean_lengths,
     disp_lengths,
     lambdas,
     lambda_median,
     lambda_OU_median,
     lambdas_OU,
     lambdas_OU_sim_low,
     lambdas_OU_sim_high,
     lambdas_sim_low,
     lambdas_sim_high,
     alphas,
     alphas_OU_sim_low,
     alphas_OU_sim_high,
     alpha_median,
     n.diffexp,
     tree_original,
     tree_3rep,
     file = file.path(simus_directory, "simulation_parameters.RData"))

