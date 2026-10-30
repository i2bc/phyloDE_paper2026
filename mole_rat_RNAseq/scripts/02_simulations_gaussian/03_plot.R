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

# Program to compute scores of the differential expression analysis methods.
# Use input from script 04_simulations_data_analysis_results.R
# The full result table from our own run is provided for convenience

## DESCRIPTION #################################################################
# Program to plot the data
################################################################################

## INPUT #######################################################################
# Formatted result table
# Produced by `mole_rat_RNAseq/scripts/02_simulations_gaussian/02_analysis.R`
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_YYY-MM-DD_gaussian_simulations_results`
# The full result table from our own run is provided for convenience
################################################################################

## OUTPUT ######################################################################
# Figure 4
################################################################################

################################################################################
################################################################################

library(ape)
library(ggplot2)
library(dplyr)
library(here)
library(cowplot)

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

datestamp_day <- "2026-01-15"
results_directory <- paste0(datestamp_day, "_", datestamp_day_simu, "_gaussian_simulations_results")
results_directory <- file.path(here_dir, "mole_rat_RNAseq/results", results_directory)
result_file <- file.path(results_directory, "result_table.tsv")

################################################################################
## Parameters
################################################################################
load(file.path(simus_directory, "simulation_parameters.RData"))

################################################################################
## Tree
################################################################################
tree <- read.tree(file = file.path(here_dir, "mole_rat_RNAseq/data/Daunesse2025/species_tree_mya.nwk"))
tree$edge.length <- tree$edge.length / max(diag(vcv(tree)))

################################################################################
## Utility function
################################################################################
get_lambda_error_BMOU <- function(phy, sigma2, sigma2_error, alpha) {
  if (alpha > 0) {
    tree_model <- phylolm::transf.branch.lengths(phy, "OUfixedRoot",
                                                 parameters = list(alpha = alpha))$tree
    tilde_t <- phyloDE:::tree_height(tree_model) / (2 * alpha)
  } else {
    tilde_t <- phyloDE:::tree_height(phy)
  }
  lambda_error <- phyloDE:::get_lambda_error(sigma2, sigma2_error, tilde_t)
  return(lambda_error)
}

###############################################################################
## Plot tree
################################################################################
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
## Colors
################################################################################
library(viridisLite)
color_methods <- mako(4, alpha = 1, begin = 0.1, end = 0.8, direction = 1)[1:4]
color_methods[4] <- inferno(4, alpha = 1, begin = 0, end = 0.8, direction = 1)[4]


################################################################################
## Load and format result table for plotting
################################################################################
## Result table
result.table <- read.table(result_file, row.names = NULL)

result.table$NP <- result.table$TP + result.table$FP

result.table[["time (log s)"]] <- log(result.table$time)
# plot_scores <- c("MCC", "FDR", "TPR", "time (log s)")#, "FP", "TP")
plot_scores <- c("FDR", "TPR")
# plot_scores <- c("MCC", "FDR", "TPR", "FP", "TP")

result.table_plot <- reshape2::melt(result.table,
                                    measure.vars = plot_scores,
                                    variable.name = "score",
                                    value.name = "score_value")
result.table_plot$method <- factor(result.table_plot$method, levels = c("limma", "eve", "phylolmBM", "phylolmOU", "phyloDEBM", "phyloDEOU"))
result.table_plot$methodMod <- result.table_plot$method
result.table_plot$methodMod <- gsub("phylolm.*", "phylolm", result.table_plot$methodMod)
result.table_plot$methodMod <- gsub("phyloDE.*", "phyloDE", result.table_plot$methodMod)
result.table_plot$methodMod <- factor(result.table_plot$methodMod, levels = c("limma", "eve", "phylolm", "phyloDE"))
result.table_plot$fitprocess <- "OU" 
result.table_plot$fitprocess[grepl("BM", result.table_plot$method)] <- "BM"
result.table_plot$fitprocess <- factor(result.table_plot$fitprocess, levels = c("OU", "BM"))

result.table_plot$Design <- result.table_plot$branch
result.table_plot$threshold <- NA
result.table_plot$threshold[result.table_plot$score == "FDR"] <- 0.05
result.table_plot$threshold[result.table_plot$score == "DR"] <- 1
result.table_plot$threshold[result.table_plot$score == "logDR"] <- 0

columnwidth <- 3 # in inches # 1 point / 72 = 1 inch

for (snr_plot in all_snr) {
  color_palete <- scale_fill_manual(values = color_methods, guide = "none", name = "")
  fill_palete <- scale_color_manual(values = color_methods, guide = "none", name = "")
  
  ###############################################################################
  ## Figure 5: OU no BM lambda and alpha three rep
  ################################################################################
  snr_fixed <- snr_plot
  lambda_error_sim_fixed <- lambda_error_base
  
  df <- subset(result.table_plot,
               snr == snr_fixed
               & tree_type == "three_rep"
               & score %in% plot_scores
               & fitprocess != "BM"
               # & Design == "MR"
               # & lambda_error_sim == lambda_error_sim_fixed
  )
  # df <- df[!(grepl("BM", df$method)), ]
  
  if (anyNA(df$score_value)) message("Some NAs in the data")
  dodge <- position_dodge(width = 0)

  df1 <- subset(df, lambda_error_sim == lambda_error_sim_fixed)
  p1 <- ggplot(df1, aes(x = log(alpha_sim), y = score_value, fill = methodMod, color = methodMod)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5,# show.legend = FALSE,
                 position = "identity", width = 0.2,
                 aes(group = interaction(method, alpha_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    color_palete + fill_palete +
    xlab("log(alpha)") +
    ylab("") +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  
  df2 <- df[sapply(df$alpha_sim, function(x) isTRUE(all.equal(x, alpha_base))), ]
  p2 <- ggplot(df2, aes(x = lambda_error_sim, y = score_value, fill = methodMod, color = methodMod)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5, #show.legend = FALSE,
                 position = "identity", #width = 0.2,
                 aes(group = interaction(method, lambda_error_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    color_palete + fill_palete +
    xlab("lambda") +
    ylab("") +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  
  prow <- plot_grid(p1, p2, labels = c("A", "B"))
  legend_b <- get_legend(
    p1 + 
      guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
  
  # colorblindr::cvd_grid()
  ggsave(file.path(results_directory, paste0(datestamp_day, "_OU_noBM_three_rep_snr_", snr_fixed, "_lambda_and_alpha.pdf")),
         height = columnwidth,
         width = 2*columnwidth,
         unit = "in")
  
  ## with tree
  top_row <- plot_grid(ggdraw(ptreeHgla), ggdraw(ptreeMR))
  pbis1 <- p1 + theme(plot.margin = unit(c(-0.3,0.1,0.1,-0.5), "cm"))
  pbis2 <- p2 + theme(plot.margin = unit(c(-0.3,0.1,0.1,-0.5), "cm"))
  pbis1 <- plot_grid(top_row, pbis1, ncol = 1, rel_heights = c(1, 6))
  pbis2 <- plot_grid(top_row, pbis2, ncol = 1, rel_heights = c(1, 6))
  
  prow <- plot_grid(pbis1, pbis2, labels = c("A", "B"))
  legend_b <- get_legend(
    p1 + 
      guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
  
  ggsave(file.path(results_directory, paste0(datestamp_day, "_with_tree_OU_noBM_three_rep_snr_", snr_fixed, "_lambda_and_alpha.pdf")),
         height = 1.3*columnwidth,
         width = 2*columnwidth,
         unit = "in")
  
  ###############################################################################
  ## extra figure: OU lambda fixed
  ################################################################################
  snr_fixed <- snr_plot
  lambda_error_sim_fixed <- lambda_error_base
  
  df <- subset(result.table_plot,
               snr == snr_fixed
               # & tree_type == "original"
               & score %in% plot_scores
               # & Design == "MR"
               & lambda_error_sim == lambda_error_sim_fixed
  )
  # df <- df[!(grepl("BM", df$method)), ]
  
  if (anyNA(df$score_value)) message("Some NAs in the data")
  
  dodge <- position_dodge(width = 0)
  
  ggplot(df, aes(x = log(alpha_sim), y = score_value, fill = methodMod, color = methodMod, linetype = fitprocess)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design, tree_type), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
                 position = "identity", width = 0.2,
                 aes(group = interaction(method, alpha_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
    scale_colour_viridis_d(end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
    xlab("log(alpha)") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          # legend.position = c(0.01, 0.01),
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  # colorblindr::cvd_grid()
  ggsave(file.path(results_directory, paste0(datestamp_day, "_OU_snr_", snr_fixed, "_lambda_", lambda_error_sim_fixed, ".pdf")),
         height = columnwidth,
         width = 2 * columnwidth,
         unit = "in")
  
  ###############################################################################
  ## extra figure: OU lambda fixed - noBM
  ################################################################################
  snr_fixed <- snr_plot
  lambda_error_sim_fixed <- lambda_error_base
  
  df <- subset(result.table_plot,
               snr == snr_fixed
               # & tree_type == "original"
               & score %in% plot_scores
               # & Design == "MR"
               & fitprocess == "OU"
               & lambda_error_sim == lambda_error_sim_fixed
  )
  # df <- df[!(grepl("BM", df$method)), ]
  
  if (anyNA(df$score_value)) message("Some NAs in the data")
  
  dodge <- position_dodge(width = 0)
  
  ggplot(df, aes(x = log(alpha_sim), y = score_value, fill = methodMod, color = methodMod, linetype = tree_type)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
                 position = "identity", width = 0.2,
                 aes(group = interaction(method, alpha_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
    scale_colour_viridis_d(end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
    xlab("log(alpha)") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          # legend.position = c(0.01, 0.01),
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  # colorblindr::cvd_grid()
  ggsave(file.path(results_directory, paste0(datestamp_day, "_OU_noBM_snr_", snr_fixed, "_lambda_", lambda_error_sim_fixed, ".pdf")),
         height = columnwidth,
         width = 2 * columnwidth,
         unit = "in")
  
  ###############################################################################
  ## extra figure: OU - alpha fixed
  ################################################################################
  snr_fixed <- snr_plot
  alpha_sim_sim_fixed <- alpha_base
  
  df <- subset(result.table_plot,
               snr == snr_fixed
               # & tree_type == "original"
               & score %in% plot_scores
               # & Design == "MR"
               & sapply(result.table_plot$alpha_sim, function(x) isTRUE(all.equal(x, alpha_base)))
  )
  
  if (anyNA(df$score_value)) message("Some NAs in the data")
  
  dodge <- position_dodge(width = 0)
  
  ggplot(df, aes(x = lambda_error_sim, y = score_value, fill = methodMod, color = methodMod, linetype = fitprocess)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design, tree_type), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
                 position = "identity", #width = 0.2,
                 aes(group = interaction(method, lambda_error_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
    scale_colour_viridis_d(end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
    xlab("lambda") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          # legend.position = c(0.01, 0.01),
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  # colorblindr::cvd_grid()
  ggsave(file.path(results_directory, paste0(datestamp_day, "_OU_snr_", snr_fixed, "_alpha_", alpha_sim_sim_fixed, ".pdf")),
         height = columnwidth,
         width = 2 * columnwidth,
         unit = "in")
  
  ###############################################################################
  ## extra figure: OU - alpha fixed - noBM
  ################################################################################
  snr_fixed <- snr_plot
  alpha_sim_sim_fixed <- alpha_base
  
  df <- subset(result.table_plot,
               snr == snr_fixed
               # & tree_type == "original"
               & score %in% plot_scores
               # & Design == "MR"
               & fitprocess == "OU"
               & sapply(result.table_plot$alpha_sim, function(x) isTRUE(all.equal(x, alpha_base)))
  )
  
  if (anyNA(df$score_value)) message("Some NAs in the data")
  
  dodge <- position_dodge(width = 0)
  
  ggplot(df, aes(x = lambda_error_sim, y = score_value, fill = methodMod, color = methodMod, linetype = tree_type)) +
    geom_hline(aes(yintercept = threshold), linetype = 2) +
    facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = FALSE,
                 position = "identity", #width = 0.2,
                 aes(group = interaction(method, lambda_error_sim, score, tree_type, Design))) +
    stat_summary(fun = median, geom = 'line', aes(group = interaction(method, score, tree_type, Design))) + 
    scale_fill_viridis_d(guide = "none", end = 0.8, option = "B") +
    scale_colour_viridis_d(end = 0.8, option = "B", guide = guide_legend(override.aes = list(size = 2))) +
    xlab("lambda") +
    ylab("") +
    theme_bw() +
    theme(text = element_text(size = 9),
          title = element_text(size = 7),
          strip.text.y = element_text(vjust = 0),
          strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          # legend.position = c(0.01, 0.01),
          # legend.justification = c("left", "bottom"),
          legend.key.size = unit(10, 'pt'),
          plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  # colorblindr::cvd_grid()
  ggsave(file.path(results_directory, paste0(datestamp_day, "_OU_noBM_snr_", snr_fixed, "_alpha_", alpha_sim_sim_fixed, ".pdf")),
         height = columnwidth,
         width = 2 * columnwidth,
         unit = "in")
}
