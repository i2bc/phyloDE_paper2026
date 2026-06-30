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
# Program to calibrate trim values
################################################################################

## INPUT #######################################################################
# Formatted tree
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Figure S9
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_trim_simulations`
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
datestamp_day <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day, "_trim_simulations")
results_directory <- file.path(here(), "mole_rat_RNAseq/results", results_directory)
dir.create(results_directory, showWarnings = FALSE)
result_file <- file.path(results_directory, "result_table.tsv")

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
nrep <- 10
ngenes <- 10000
all_branch <- c("MR", "Hgla")

all_lambda_error <- c(0.1, 0.25, 0.5, 0.75, 0.9)
lambda_error_base <- 0.75
all_alpha <- log(2) / c(0.05, 0.1, 0.5, 1, 5)
alpha_base <- log(2) / 0.5

## params OU
sapply(all_alpha, function(x) sapply(all_lambda_error, function(y) get_parameters_OU(x, y, tree_original)))

################################################################################
## Trim
################################################################################

all_trim_test <- list(c(0.30, 0.00),
                      c(0.25, 0.05),
                      c(0.20, 0.10),
                      c(0.15, 0.15))

################################################################################
## Simulations and analysis
################################################################################
all_results <- NULL

for (tree_type in all_trees) {
  tree_sim <- switch(tree_type, original = tree_original, three_rep = tree_3rep)
  colDataTree <- colData[match(tree_sim$tip.label, rownames(colData)), ]
  for (branch in all_branch) {
    for (lambda_error_sim in all_lambda_error) {
      for (alpha_sim in all_alpha) {
        if ((lambda_error_sim != lambda_error_base) && (alpha_sim != alpha_base)) next
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
          ## Normalize
          total_variance <- get_total_variance_OU(phy = tree_sim, sigma2 = params_sim$sigma2_tree, sigma2_error = params_sim$sigma2_error, alpha = alpha_sim)
          sim_data <- sim_data / sqrt(total_variance)
          
          ################################################################################
          ## phyloDE
          ################################################################################
          ## Design
          if (branch == "Hgla") {
            condition <- colDataTree$conditionHgla
          } else if (branch == "MR") {
            condition <- colDataTree$conditionMR
          }
          design_formula <- as.formula(~ condition)
          design_data <- data.frame(condition = condition)
          rownames(design_data) <- colnames(sim_data)
          design_data$condition <- factor(design_data$condition)
          design <- model.matrix(design_formula, design_data)
          
          ################################################################################
          ## phyloDE - BM
          fitBM <- phyloDE::phylolmFit(sim_data,
                                       design = design,
                                       phy = tree_sim,
                                       model = "BM",
                                       measurement_error = TRUE,
                                       use_consensus = TRUE,
                                       ddf_method = "Samples",
                                       ncores = 6)
          
          ## lambda 
          lambdasBM <- tanh(fitBM$consensus_tree$params$atanh_lambda)
          
          ################################################################################
          ## phyloDE BM - all trim values
          for (trim_lambda in all_trim_test) {
            rho <- tanh(phyloDE:::mean_trim(atanh(lambdasBM), trim = trim_lambda))
            consensus_tree <- fitBM$consensus_tree
            consensus_tree$tree <- phylolm::transf.branch.lengths(tree_sim, model = "lambda",
                                                                  parameters = list(lambda = rho))$tree
            consensus_tree$params$tree_ind <- NULL
            
            fitBM <- phyloDE::phylolmFit(sim_data,
                                         design = design,
                                         phy = tree_sim,
                                         model = 'BM',
                                         measurement_error = TRUE,
                                         use_consensus = TRUE,
                                         consensus_tree = consensus_tree,
                                         ddf_method = "Samples")
            fitBM <- limma::eBayes(fitBM, trend = FALSE)
            BMPvalues <- fitBM$p.value[, ncol(fitBM$p.value)]
            
            resultBM <- data.frame(branch = branch,
                                   tree_type = tree_type,
                                   alpha = alpha_sim,
                                   sigma2_error = params_sim$sigma2_error,
                                   sigma2_tree = params_sim$sigma2_tree,
                                   lambda = lambda_error_sim,
                                   typeI = mean(BMPvalues <= 0.05),
                                   typeIadj = mean(p.adjust(BMPvalues, method = 'BH') <= 0.05),
                                   method = "BM",
                                   trim = paste0(trim_lambda, collapse = "_"),
                                   trim_alpha = NA)
            
            all_results <- rbind(all_results, resultBM)
            write.table(all_results,
                        result_file,
                        sep='\t', quote = FALSE)
          }
          
          ################################################################################
          ## phyloDE - OU
          fitOU <- phyloDE::phylolmFit(sim_data,
                                       design = design,
                                       phy = tree_sim,
                                       model = "OUfixedRoot",
                                       measurement_error = TRUE,
                                       use_consensus = TRUE,
                                       ddf_method = "Samples",
                                       ncores = 6)
          
          ## alpha and lambda
          rho_atanh <- fitOU$consensus_tree$params$trans_alpha
          alpha_log <- fitOU$consensus_tree$params$log_alpha
          lambda_OU_atanh <- fitOU$consensus_tree$params$atanh_lambda_error
          
          ################################################################################
          ## phyloDE - OU - rho
          for (trim_lambda in all_trim_test) {
            for (trim_alpha in all_trim_test) {
              trans_inv_alpha <- function(tt) {
                phyloDE:::rho_prime_inv(1 - tanh(tt), t_original_tree, c(0, 300))
              }
              alpha_bis <- trans_inv_alpha(phyloDE:::mean_trim(rho_atanh, trim = trim_alpha))
              lambda_bis <- tanh(phyloDE:::mean_trim(lambda_OU_atanh, trim = trim_lambda, na.rm = TRUE))
              
              tree_model <- phylolm::transf.branch.lengths(tree_sim, "OUfixedRoot", parameters = list(alpha = alpha_bis))$tree
              tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_bis))$tree
              tree_model <- phyloDE:::rescale_tree(tree_model)
              consensus_tree <- fitOU$consensus_tree
              consensus_tree$tree <- tree_model
              consensus_tree$params$tree_ind <- NULL
              
              fitOU <- phyloDE::phylolmFit(sim_data,
                                           design = design,
                                           phy = tree_sim,
                                           model = 'OUfixedRoot',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE, 
                                           consensus_tree = consensus_tree,
                                           ddf_method = "Samples")
              fitOU <- limma::eBayes(fitOU, trend = FALSE)
              OUPvalues <- fitOU$p.value[, ncol(fitOU$p.value)]
              
              resultOU <- data.frame(branch = branch,
                                     tree_type = tree_type,
                                     alpha = alpha_sim,
                                     sigma2_error = params_sim$sigma2_error,
                                     sigma2_tree = params_sim$sigma2_tree,
                                     lambda = lambda_error_sim,
                                     typeI = mean(OUPvalues <= 0.05),
                                     typeIadj = mean(p.adjust(OUPvalues, method = 'BH') <= 0.05),
                                     method = "OU",
                                     trim = paste0(trim_lambda, collapse = "_"),
                                     trim_alpha = paste0(trim_alpha, collapse = "_"))
              
              all_results <- rbind(all_results, resultOU)
              write.table(all_results,
                          result_file,
                          sep='\t', quote = FALSE)
            }
          }
          
          ################################################################################
          ## phyloDE - OU - alpha
          for (trim_lambda in all_trim_test) {
            for (trim_alpha in all_trim_test) {
              alpha_bis <- exp(phyloDE:::mean_trim(alpha_log, trim = trim_alpha))
              lambda_bis <- tanh(phyloDE:::mean_trim(lambda_OU_atanh, trim = trim_lambda, na.rm = TRUE))
              
              tree_model <- phylolm::transf.branch.lengths(tree_sim, "OUfixedRoot", parameters = list(alpha = alpha_bis))$tree
              tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda_bis))$tree
              tree_model <- phyloDE:::rescale_tree(tree_model)
              consensus_tree <- fitOU$consensus_tree
              consensus_tree$tree <- tree_model
              consensus_tree$params$tree_ind <- NULL
              
              fitOU <- phyloDE::phylolmFit(sim_data,
                                           design = design,
                                           phy = tree_sim,
                                           model = 'OUfixedRoot',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE, 
                                           consensus_tree = consensus_tree,
                                           ddf_method = "Samples")
              fitOU <- limma::eBayes(fitOU, trend = FALSE)
              OUPvalues <- fitOU$p.value[, ncol(fitOU$p.value)]
              
              resultOU <- data.frame(branch = branch,
                                     tree_type = tree_type,
                                     alpha = alpha_sim,
                                     sigma2_error = params_sim$sigma2_error,
                                     sigma2_tree = params_sim$sigma2_tree,
                                     lambda = lambda_error_sim,
                                     typeI = mean(OUPvalues <= 0.05),
                                     typeIadj = mean(p.adjust(OUPvalues, method = 'BH') <= 0.05),
                                     method = "OU_alpha",
                                     trim = paste0(trim_lambda, collapse = "_"),
                                     trim_alpha = paste0(trim_alpha, collapse = "_"))
              
              all_results <- rbind(all_results, resultOU)
              write.table(all_results,
                          result_file,
                          sep='\t', quote = FALSE)
            }
          }
        }
      }
    }
  }
}

################################################################################
################################################################################
## plot
################################################################################
################################################################################

datestamp_day <- "2026-02-17"
data_file <- "mole_rat_RNAseq/data"
results_directory <- paste0(datestamp_day, "_trim_simulations")
results_directory <- file.path(here(), "mole_rat_RNAseq/results", results_directory)
result_file <- file.path(results_directory, "result_table.tsv")
all_results <- read.table(result_file)

################################################################################
## colors
library(ggplot2)
library(cowplot)
library(viridisLite)
color_methods <- mako(4, alpha = 1, begin = 0.1, end = 0.8, direction = 1)[1:4]
color_methods[4] <- inferno(4, alpha = 1, begin = 0, end = 0.8, direction = 1)[4]
color_palete <- scale_fill_manual(values = color_methods[c(1, 3, 4)], guide = "none")
fill_palete <- scale_color_manual(values = color_methods[c(1, 3, 4)], guide = "none")
columnwidth <- 3 # in inches # 1 point / 72 = 1 inch

###############################################################################
## Plot tree
# load tree and drop replicates
tree <- tree_3rep
tree <- drop.tip(tree, tree$tip.label[duplicated(sub("_.*$", "", tree$tip.label))])
tree$tip.label <- sub("_.*$", "", tree$tip.label)

col_orig <- c("grey60", "black")

ptreeHgla <- function() {
  par(mar = c(0, 2.5, 0, 0))
  col_plot <-  rep(col_orig[1], 6)
  col_plot[3] <- col_orig[2]
  col_tip <-  rep(col_orig[1], 4)
  col_tip[3] <- col_orig[2]
  plot(tree, edge.color = col_plot, edge.width = 3,
       y.lim = c(-0.1, 4.5), x.lim = c(0, 1.46),
       show.tip.label = TRUE, cex = 0.6, label.offset = 0.1,
       tip.color = col_tip)
  #
}
# ggdraw(ptreeHgla)

ptreeMR <- function() {
  par(mar = c(0, 1.6, 0, 0.9))
  col_plot <-  rep(col_orig[1], 6)
  col_plot[c(2, 3, 4)] <- col_orig[2]
  col_tip <-  rep(col_orig[1], 4)
  col_tip[c(2, 3)] <- col_orig[2]
  plot(tree, edge.color = col_plot, edge.width = 3,
       y.lim = c(-0.1, 4.5), x.lim = c(0, 1.46),
       show.tip.label = TRUE, cex = 0.6, label.offset = 0.1,
       tip.color = col_tip)
}

###############################################################################
## Plot figure S8

all_results <- subset(all_results, method == "BM" | trim == trim_alpha)

all_results$Design <- all_results$branch
all_results$trim <- sub("_", ", ", all_results$trim)
all_results <- subset(all_results, !(trim == "0.3, 0") & tree_type == "three_rep")
all_results$method <- factor(all_results$method)
levels(all_results$method) <- c("BM", "OU (rho)", "OU (alpha)")

df1 <- subset(all_results, lambda == lambda_error_base)
p1 <- ggplot(df1, aes(x = log(alpha), y = typeI, fill = trim, color = trim, group = interaction(trim, method, lambda, branch))) +
  geom_hline(aes(yintercept = 0.05), linetype = 2) +
  facet_grid(cols = vars(Design), rows = vars(method), scales = "free", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5,# show.legend = FALSE,
               position = "identity",
               aes(group = interaction(alpha, trim, method, lambda, branch))) +
  stat_summary(fun = median, geom = 'line') + 
  color_palete + fill_palete +
  xlab("log(alpha)") +
  ylab("type I error") +
  ylim(c(0, 0.2)) +
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

df2 <- all_results[sapply(all_results$alpha, function(x) isTRUE(all.equal(x, alpha_base))), ]
p2 <- ggplot(df2, aes(x = lambda, y = typeI, fill = trim, color = trim, group = interaction(trim, method, alpha, branch))) +
  geom_hline(aes(yintercept = 0.05), linetype = 2) +
  facet_grid(cols = vars(Design), rows = vars(method), scales = "free", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5,# show.legend = FALSE,
               position = "identity",
               aes(group = interaction(alpha, trim, method, lambda, branch))) +
  stat_summary(fun = median, geom = 'line') + 
  color_palete + fill_palete +
  xlab("lambda") +
  ylab("type I error") +
  ylim(c(0, 0.2)) +
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

prow <- plot_grid(p1, p2, labels = c("A", "B"))
legend_b <- get_legend(
  p1 + 
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

ggsave(file.path(results_directory, "all_plot.pdf"),
       height = 2*columnwidth,
       width = 2*columnwidth,
       unit = "in")

## with tree
top_row <- plot_grid(ggdraw(ptreeHgla), ggdraw(ptreeMR))
pbis1 <- p1 + theme(plot.margin = unit(c(-0.3,0.1,0.1,0.1), "cm"))
pbis2 <- p2 + theme(plot.margin = unit(c(-0.3,0.1,0.1,0.1), "cm"))
pbis1 <- plot_grid(top_row, pbis1, ncol = 1, rel_heights = c(1, 6))
pbis2 <- plot_grid(top_row, pbis2, ncol = 1, rel_heights = c(1, 6))

prow <- plot_grid(pbis1, pbis2, labels = c("A", "B"))
legend_b <- get_legend(
  p1 + 
    guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

ggsave(file.path(results_directory, "all_plot_with_tree.pdf"),
       height = 2*columnwidth,
       width = 2*columnwidth,
       unit = "in")
