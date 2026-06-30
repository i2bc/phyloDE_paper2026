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
# Program to analyse a single simulation on 10000 genes
################################################################################

## INPUT #######################################################################
# Formatted tree
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Figures 1, S6, S7, S8
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

tree_sim <- tree_3rep

## tree without replicates
tree <- read.tree(file = here("mole_rat_RNAseq/data/data_formatted/species_tree_mya.nwk"))
tree$edge.length <- tree$edge.length / max(diag(vcv(tree)))

################################################################################
## Result files
################################################################################
here_dir <- here()
datestamp_day <- "2026-03-04" #format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day, "_test_one_simulation")
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

get_sigma2_error <- function(par_eve) {
  par_eve[,"sigma2"] / 2 / par_eve[,"alpha"] * par_eve[,"beta"]
}

get_rho <- function(alpha, t_tree) {
  1 - (1 - exp(-2 * alpha * t_tree)) / (2 * alpha * t_tree)
}

################################################################################
## Parameters
################################################################################
nrep <- 10
ngenes <- 10000
all_branch <- c("MR", "Hgla")

lambda_error_sim <- 0.75
alpha_sim <- log(2) / 0.5
rho_sim <- get_rho(alpha_sim, t_original_tree)

################################################################################
################################################################################
## Simulations with fixed parameters - no size effect - mean 1
################################################################################
################################################################################

params_sim <- get_parameters_OU(alpha_sim, lambda_error_sim, tree_sim)
beta_sim <- 2 * alpha_sim * params_sim$sigma2_error / params_sim$sigma2_tree
set.seed(1289)
################################################################################
## Simulate OU with fixed parameters
sim_data <- t(phylolm::rTrait(ngenes, phy = tree_sim,
                              model = "OU", parameters = list(sigma2 = params_sim$sigma2_tree, alpha = alpha_sim, ancestral.state = 1, optimal.value = 1)))
rownames(sim_data) <- paste0("g", 1:nrow(sim_data))
################################################################################
## Add noise
sim_data <- sim_data + matrix(rnorm(prod(dim(sim_data)), 0, sqrt(params_sim$sigma2_error)), ncol = ncol(sim_data))
################################################################################
## Normalize
total_variance <- get_total_variance_OU(phy = tree_sim, sigma2 = params_sim$sigma2_tree, sigma2_error = params_sim$sigma2_error, alpha = alpha_sim)
sim_data <- sim_data / sqrt(total_variance)

################################################################################
################################################################################
## Analysis with phyloDE and eve - no size effect
################################################################################
################################################################################

################################################################################
## phyloDE - BM
fitBM <- phyloDE::phylolmFit(sim_data,
                             design = NULL,
                             phy = tree,
                             model = "BM",
                             measurement_error = TRUE,
                             use_consensus = TRUE,
                             ddf_method = "Samples",
                             ncores = 6)

## lambda 
lambda_BM <- tanh(fitBM$consensus_tree$params$atanh_lambda)
hist(lambda_BM)

res <- data.frame(parameters = "lambda",
                  process = "BM", 
                  value = lambda_BM,
                  true = lambda_error_sim,
                  trim_mean = fitBM$consensus_tree$params$lambda_error,
                  type = "transformed")

################################################################################
## phyloDE - OU
fitOU <- phyloDE::phylolmFit(sim_data,
                             phy = tree,
                             model = "OUfixedRoot",
                             measurement_error = TRUE,
                             use_consensus = TRUE,
                             ddf_method = "Samples",
                             ncores = 6)

## alpha and lambda
rho_OU <- tanh(fitOU$consensus_tree$params$trans_alpha)
alpha_OU <- exp(fitOU$consensus_tree$params$log_alpha)
lambda_OU <- tanh(fitOU$consensus_tree$params$atanh_lambda_error)
sigma2_OU <- fitOU$sigma^2
sigma_err_OU <- sigma2_OU * (1 - lambda_OU) / lambda_OU
estim_alpha_OU <- fitOU$consensus_tree$params$alpha
estim_rho_OU <- 1 - (1 - exp(-2 * estim_alpha_OU * t_original_tree)) / (2 * estim_alpha_OU * t_original_tree)
hist(rho_OU)
hist(alpha_OU)
hist(lambda_OU)
theta_OU <- fitOU$coefficients[, 1]

## get log likelihood
ll_OU <- rep(NA, nrow(sim_data))
for (i in 1:nrow(sim_data)) {
  dd <- sim_data[i, ]
  names(dd) <- colnames(sim_data)
  tree_bis <- phylolm::transf.branch.lengths(tree_sim, model = "OUfixedRoot", parameters = list(alpha = alpha_OU[i]))$tree
  ll_OU[i] <- phylolm::phylolm(dd ~ 1, phy = tree_bis, model = "lambda", measurement_error = FALSE,
                               lower.bound = list(lambda = lambda_OU[i]),
                               upper.bound = list(lambda = lambda_OU[i]),
                               starting.value = list(lambda = lambda_OU[i]))$logLik
}

res <- rbind(res,
             data.frame(parameters = "rho",
                        process = "OU", 
                        value = rho_OU,
                        true = rho_sim,
                        trim_mean = estim_rho_OU,
                        type = "transformed"),
             data.frame(parameters = "lambda",
                        process = "OU", 
                        value = lambda_OU,
                        true = lambda_error_sim,
                        trim_mean = fitOU$consensus_tree$params$lambda_error,
                        type = "transformed"),
             data.frame(parameters = "log(alpha)",
                        process = "OU", 
                        value = log(alpha_OU),
                        true = log(alpha_sim),
                        trim_mean = NA,
                        type = "natural"),
             data.frame(parameters = "log(s^2)",
                        process = "OU", 
                        value = log(sigma_err_OU),
                        true = log(params_sim$sigma2_error),
                        trim_mean = NA,
                        type = "natural")
             )

################################################################################
## EVE
colSpecies <- sub("_.*$","", colnames(sim_data))

fitEVE <- evemodel::fitOneTheta(tree = tree, gene.data = sim_data,
                                colSpecies = colSpecies,
                                cores = 6)
## alpha and lambda
alpha_eve <- fitEVE$par[, "alpha"]
rho_eve <- get_rho(fitEVE$par[, "alpha"], 1)
sigma2_err_eve <- get_sigma2_error(fitEVE$par)
lambda_eve <- phyloDE:::get_lambda_error(fitEVE$par[, "sigma2"], sigma2_err_eve, 1)
beta_eve <- fitEVE$par[, "beta"]
theta_eve <- fitEVE$par[, "theta"]
hist(rho_eve)
hist(alpha_eve)
hist(lambda_eve)
ll_EVE <- fitEVE$ll

res <- rbind(res,
             data.frame(parameters = "rho",
                        process = "EVE", 
                        value = rho_eve,
                        true = rho_sim,
                        trim_mean = NA,
                        type = "transformed"),
             data.frame(parameters = "lambda",
                        process = "EVE", 
                        value = lambda_eve,
                        true = lambda_error_sim,
                        trim_mean = NA,
                        type = "transformed"),
             data.frame(parameters = "log(alpha)",
                        process = "EVE", 
                        value = log(alpha_eve),
                        true = log(alpha_sim),
                        trim_mean = NA,
                        type = "natural"),
             data.frame(parameters = "log(beta)",
                        process = "EVE", 
                        value = log(fitEVE$par[, "beta"]),
                        true = log(beta_sim),
                        trim_mean = NA,
                        type = "natural")
)

## compare parameters and likelihoods
all.equal(ll_OU, ll_EVE)
get_rho(range(alpha_OU), 1)
get_rho(range(alpha_eve), 1)

range(fitEVE$par[, "beta"])
range(fitEVE$par[, "alpha"])

plot(ll_OU, ll_EVE)
plot(rho_OU, rho_eve)
plot(lambda_OU, lambda_eve)
plot(theta_OU, theta_eve)

hist(theta_OU - theta_eve)
hist(ll_OU - ll_EVE)

################################################################################
################################################################################
## plots
################################################################################
################################################################################

## colors
library(ggplot2)
library(cowplot)
library(viridisLite)
color_methods <- mako(4, alpha = 1, begin = 0.1, end = 0.8, direction = 1)[1:4]
color_methods[4] <- inferno(4, alpha = 1, begin = 0, end = 0.8, direction = 1)[4]
columnwidth <- 3 # in inches # 1 point / 72 = 1 inch


################################################################################
## Figure S6 - phyloDE OU

resOU <- subset(res, process == "OU")
resOU$method <- "phylolm"
resOU$parameters <- factor(resOU$parameters, levels = c("log(s^2)", "log(alpha)", "lambda", "rho"))
pOU <- ggplot(resOU, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.5, color = color_methods[1], fill = color_methods[1]) +
  facet_wrap(vars(parameters), nrow = 1, scale = "free", strip.position = "bottom", labeller = label_parsed) +
  # ggh4x::facet_grid2(rows = vars(method), cols = vars(parameters), scale = "free", switch = "x", independent = "y") +
  geom_vline(aes(xintercept = true), color = color_methods[3], linewidth = 1) +
  geom_vline(aes(xintercept = trim_mean), linetype = "dashed", color = color_methods[4], linewidth = 0.9) +
  labs(x = NULL) +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        # legend.justification = c("left", "bottom"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
pOU

ggsave(file.path(results_directory, "OU_param_estim_plot.pdf"),
       height = 0.6*columnwidth,
       width = 2*columnwidth,
       unit = "in")

################################################################################
## figure S6 - EVE

resEVE <- subset(res, process == "EVE")
resEVE$method <- "eve"
resEVE$parameters <- factor(resEVE$parameters, levels = c("log(beta)", "log(alpha)", "lambda", "rho"))
pEVE <- ggplot(resEVE, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, alpha = 0.5, color = color_methods[1], fill = color_methods[1]) +
  facet_wrap(vars(parameters), nrow = 1, scale = "free", strip.position = "bottom", labeller = label_parsed) +
  # ggh4x::facet_grid2(rows = vars(method), cols = vars(parameters), scale = "free", switch = "x", independent = "y") +
  geom_vline(aes(xintercept = true), color = color_methods[3], linewidth = 1) +
  labs(x = NULL) +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        # legend.justification = c("left", "bottom"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
pEVE

ggsave(file.path(results_directory, "EVE_param_estim_plot.pdf"),
       height = 0.6*columnwidth,
       width = 2*columnwidth,
       unit = "in")

pOUEVE <- plot_grid(pOU, pEVE, ncol = 1, rel_heights = c(1, 1))
pOUEVE

ggsave(file.path(results_directory, "EVEOU_param_estim_plot.pdf"),
       height = 1.2*columnwidth,
       width = 2*columnwidth,
       unit = "in")

################################################################################
## figure 1 - EVE only naturals

resEVE$parameters <- factor(resEVE$parameters, levels = c("theta", "log(beta)", "log(alpha)"))

ggplot(resEVE, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 100, alpha = 0.5, color = color_methods[1], fill = color_methods[1]) +
  facet_wrap(vars(parameters), scale = "free", strip.position = "bottom", labeller = label_parsed, dir = "bl", ncol = 3) +
  geom_vline(aes(xintercept = true), color = color_methods[3], linewidth = 0.8) +
  labs(x = NULL) +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        # legend.justification = c("left", "bottom"),
        legend.key.size = unit(10, 'pt'),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

ggsave(file.path(results_directory, "EVE_param_natural_one_line.pdf"),
       height = 0.5*columnwidth,
       width = 1.5*columnwidth,
       unit = "in")

################################################################################
## figure S8 - bi plots
library(ggside)
library(latex2exp)

make_plot <- function(dd, titleplot) {
  pp <- ggplot(dd, aes(x = phylolm, y = eve)) + 
    geom_point(alpha = 0.1, color = color_methods[1], fill = color_methods[1]) +
    coord_equal() +
    geom_abline(slope = 1) + 
    geom_hline(aes(yintercept = truevalue), color = color_methods[3], linewidth = 0.8) +
    geom_vline(aes(xintercept = truevalue), color = color_methods[3], linewidth = 0.8) +
    geom_xsidehistogram(aes(y = after_stat(density)),
                        bins = 100, alpha = 0.5, color = color_methods[1], fill = color_methods[1]) + 
    geom_ysidehistogram(aes(x = after_stat(density)),
                        bins = 100, alpha = 0.5, color = color_methods[1], fill = color_methods[1]) + 
    geom_xsidevline(aes(xintercept = truevalue), color = color_methods[3], linewidth = 0.8) +
    geom_xsidevline(aes(xintercept = trim_mean), linetype = "dashed", color = color_methods[4], linewidth = 0.9) +
    geom_ysidehline(aes(yintercept = truevalue), color = color_methods[3], linewidth = 0.8) +
    theme_bw() +
    labs(title = titleplot) + 
    theme(text = element_text(size = 9),
          title = element_text(size = 9),
          ggside.panel.scale = .5
          # plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
    )
  # pp <- ggMarginal(pp, type="histogram",
  #                  bins = 100, alpha = 0.5, color = color_methods[1], fill = color_methods[1])
  return(pp)
}

resll <- data.frame(eve = ll_EVE, phylolm = ll_OU, truevalue = NA, trim_mean = NA)
pll <- make_plot(resll, "Optimized log likelihood")
pll

restheta <- data.frame(eve = theta_eve, phylolm = theta_OU, truevalue = 1, trim_mean = NA)
ptheta <- make_plot(restheta, TeX(r"(Parameter $\theta$)"))
ptheta

resrho <- data.frame(eve = rho_eve, phylolm = rho_OU, truevalue = rho_sim, trim_mean = estim_rho_OU)
prho <- make_plot(resrho, TeX(r"(Parameter $\rho$)"))
prho <- prho + coord_equal(xlim = c(0, 1), ylim = c(0, 1))
prho

reslambda <- data.frame(eve = lambda_eve, phylolm = lambda_OU, truevalue = lambda_error_sim, trim_mean = fitOU$consensus_tree$params$lambda_error)
plambda <- make_plot(reslambda, TeX(r"(Parameter $\lambda$)"))
plambda <- plambda + coord_equal(xlim = c(0, 1), ylim = c(0, 1))
plambda

plot_grid(pll, ptheta, prho, plambda, ncol = 2)

ggsave(file.path(results_directory, "bi_plots_eve_phyloDE.pdf"),
       height = 2*columnwidth,
       width = 2*columnwidth,
       unit = "in")
