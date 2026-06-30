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
# Program to plot the formated results
################################################################################

## INPUT #######################################################################
# Formatted result table
# Saved in `mole_rat_RNAseq/results/YYYY-MM-DD_YYYY-MM-DD_simulations_results_mole_rat_RNAseq_data_Liver`
################################################################################

## OUTPUT ######################################################################
# Analysis plots
################################################################################

################################################################################
################################################################################

library(compcodeR)
library(ggplot2)
library(cowplot)
library(dplyr)
library(here)
library(ape)

################################################################################
## File management
################################################################################
here_dir <- here()
data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"
data_name <- paste0("mole_rat_RNAseq_data_", tissue)
datestamp_day_simus <- "2026-01-19" ## Change here for the simulation date
datestamp_day_anaysis <- "2026-01-19" ## Change here for the analysis date
simus_directory <- paste0(datestamp_day_simus, "_simulations_", data_name)
simus_directory <- file.path(here_dir, "mole_rat_RNAseq/results", simus_directory)
results_directory <- paste0(datestamp_day_anaysis, "_", datestamp_day_simus, "_simulations_results_", data_name)
results_directory <- file.path(here_dir, "mole_rat_RNAseq/results", results_directory)

################################################################################
## Parameters
################################################################################
load(file.path(simus_directory, "simulation_parameters.RData"))

## Iterations
Nmin <- 1
Nmax <- 20

## Base parameters
base_effect_size <- 5
base_prop_var_tree <- 0.8
base_fact_disp <- 1
base_cond_type <- "Hgla"

################################################################################
## Load and format result table for plotting
################################################################################
## Result table
result.table <- readRDS(file.path(results_directory, paste0("full_result_table_", Nmin, "_", Nmax, ".rds")))

## Indices with no discovery : fp + tp = 0
insnodisc <- !is.na(result.table$fp) & (result.table$fp == 0) & !is.na(result.table$tp) & (result.table$tp == 0)
## Replace NA by 0 in FDR when fp + tp = 0
result.table[is.na(result.table$fdr) & insnodisc, "fdr"] <- 0.0
## Replace -1000 by 0 in FDR when no discovery
result.table[result.table$fdr == -1000 & insnodisc, "fdr"] <- 0.0
## Replace NA by 0 in MCC when fp + tp = 0
result.table[is.na(result.table$mcc) & insnodisc, "mcc"] <- 0.0
## Replace -1000 by 0 in MCC when fp + tp = 0
result.table[result.table$mcc == -1000 & insnodisc, "mcc"] <- 0.0

## Names
colnames(result.table)[colnames(result.table) == "fdr"] <- "FDR"
colnames(result.table)[colnames(result.table) == "tpr"] <- "TPR"
colnames(result.table)[colnames(result.table) == "auc"] <- "AUC"
colnames(result.table)[colnames(result.table) == "mcc"] <- "MCC"
colnames(result.table)[colnames(result.table) == "fp"] <- "FP"
colnames(result.table)[colnames(result.table) == "tp"] <- "TP"

result.table$methodMod <- sub("OUfixedRoot", "OU", result.table$methodMod)
result.table$methodMod <- sub("lengthNorm.limma", "limma", result.table$methodMod)

result.table$Design <- result.table$cond_type
# result.table$Design <- sub("rodents", "rodent", result.table$Design)
# result.table$Design <- sub("primates", "primate", result.table$Design)
#result.table$Design <- sub("bloc", "block", result.table$Design)
# result.table$Design <- factor(result.table$Design, levels = c("rodent", "primate"))

result.table$prop_var_tree <- as.factor(result.table$prop_var_tree)
result.table$effect_size <- as.factor(result.table$effect_size)
result.table$sel_strength <- as.factor(result.table$sel_strength)

## Indicator emp no emp
result.table$emp <- "fixed"
result.table$emp[result.table$prop_var_tree == "simlow"] <- "low var"
result.table$emp[result.table$prop_var_tree == "simhigh"] <- "high var"

## names
ind <- which(result.table$method == "phyloDE")
result.table[ind, ]$methodMod <- sapply(ind, function(i) paste0(result.table[i, ]$method, "_",
                                                                result.table[i, ]$methodProcess, "_",
                                                                paste0("eBayes_", result.table[i, ]$eBayes), "_",
                                                                paste0("ddf_method_", result.table[i, ]$ddf_method), "_",
                                                                result.table[i, ]$trend, "_",
                                                                result.table[i, ]$regCor))
result.table$NP <- result.table$TP + result.table$FP

## Format for plot
result.table_plot <- reshape2::melt(result.table,
                                    measure.vars = c("AUC", "MCC", "FDR", "TPR", "FP", "TP", "NP", "time"),
                                    variable.name = "score",
                                    value.name = "score_value")
result.table_plot$threshold <- NA
result.table_plot$threshold[result.table_plot$score == "FDR"] <- 0.05
result.table_plot$threshold[result.table_plot$score == "DR"] <- 1
result.table_plot$threshold[result.table_plot$score == "logDR"] <- 0
# result.table_plot$threshold[result.table_plot$score == "TD"] <- n.diffexp

pretty_names <- c(
  "limma_no_trend_cor" = "limma (notrend)",
  "limma_with_trend_cor" = "limma",
  "phylolm_BM" = "phylolm_raw (BM)",
  "phylolm_OU" = "phylolm_raw (OU)",
  "phyloDE_BM_eBayes_FALSE_ddf_method_Samples_with_trend_no_reg" = "phylolm (BM)",
  "phyloDE_OUfixedRoot_eBayes_FALSE_ddf_method_Samples_with_trend_no_reg" = "phylolm (OU)",
  "evemodel_chi2" = "evemodel",
  "evemodel_emp" = "evemodel (emp)",
  "phyloDE_BM_eBayes_TRUE_ddf_method_Samples_no_trend_reg" = "phyloDE (BM, notrend)",
  "phyloDE_OUfixedRoot_eBayes_TRUE_ddf_method_Samples_no_trend_reg" = "phyloDE (OU, notrend)",
  "phyloDE_BM_eBayes_TRUE_ddf_method_Samples_with_trend_reg" = "phyloDE (BM)",
  "phyloDE_OUfixedRoot_eBayes_TRUE_ddf_method_Samples_with_trend_reg" = "phyloDE (OU)",
  "phyloDE_BM_eBayes_FALSE_ddf_method_Samples_with_trend_reg" = "phyloDE (BM, noeBayes)",
  "phyloDE_BM_eBayes_TRUE_ddf_method_Samples_with_trend_no_reg" = "phyloDE (BM, noreg)",
  "phyloDE_OUfixedRoot_eBayes_FALSE_ddf_method_Samples_with_trend_reg" = "phyloDE (OU, noeBayes)",
  "phyloDE_OUfixedRoot_eBayes_TRUE_ddf_method_Samples_with_trend_no_reg" = "phyloDE (OU, noreg)"
)

order_methods <- names(pretty_names)

result.table_plot$methodMod <- factor(result.table_plot$methodMod, levels = order_methods)
levels(result.table_plot$methodMod) <- pretty_names[levels(result.table_plot$methodMod)]

plot_scores <- c("FDR", "TPR")#, "time")#, "FP", "TP")
plot_scores_bis <- c("FDR", "NP", "time")#, "FP", "TP")

result.table_plot$emp <- factor(result.table_plot$emp, levels = c("fixed", "low var", "high var"))

result.table_plot$methodBare <- result.table_plot$methodMod
result.table_plot$methodBare <- gsub("phylolm.*", "phylolm", result.table_plot$methodBare)
result.table_plot$methodBare <- gsub("phyloDE.*", "phyloDE", result.table_plot$methodBare)
result.table_plot$methodBare <- gsub("limma.*", "limma", result.table_plot$methodBare)
result.table_plot$methodBare <- gsub("evemodel.*", "eve", result.table_plot$methodBare)
result.table_plot$methodBare <- factor(result.table_plot$methodBare, levels = c("limma", "eve", "phylolm", "phyloDE"))
result.table_plot$fitprocess <- "OU" 
result.table_plot$fitprocess[grepl("BM", result.table_plot$methodMod)] <- "BM"
result.table_plot$fitprocess <- factor(result.table_plot$fitprocess, levels = c("OU", "BM"))

columnwidth <- 3 # in inches # 1 point / 72 = 1 inch

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
  par(mar = c(0, 1.75, 0, 0.75))
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

color_ebayes <- inferno(4, alpha = 1, begin = 0, end = 0.8, direction = 1)[c(1, 2, 4)]


###############################################################################
## Figure 2: fixed parameters effect size fixed
################################################################################
effect_size_plot <- 4
df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM")
             & use_lengths ==  "with_lengths"
             & tree_type == "three_rep"
             & transformation %in% c("log2")
             & effect_size == effect_size_plot
             & fact_disp == base_fact_disp
             & emp == "fixed"
             & score %in% plot_scores
)
df <- df[!(grepl("_raw", df$methodMod)), ]
df <- df[!(grepl("notrend", df$methodMod)), ]
df <- df[!(grepl("noeBayes", df$methodMod)), ]
df <- df[!(grepl("noreg", df$methodMod)), ]
df <- df[!is.na(df$methodMod), ]
df <- subset(df, fitprocess != "BM")

if (nrow(df) / 2 / 2 / 4 / 2  != Nmax - Nmin + 1) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

## plot data
p <- ggplot(df, aes(x = model_process, y = score_value, fill = methodBare, color = methodBare)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = TRUE,
               position = "identity", width = 0.2,
               aes(group = interaction(methodMod, score, emp, model_process, Design, tree_type))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(use_lengths, tree_type, methodMod, Design), colour = methodBare)) +
  scale_fill_manual(values = color_methods, name = "") +
  scale_color_manual(values = color_methods,  name = "") +
  xlab("Simulation Process") +
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
        legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.1,-0.5), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
p

## add tree
bottom_row <- p + theme(plot.margin = unit(c(-0.3,0.1,0.1,-0.5), "cm"))
top_row <- plot_grid(ggdraw(ptreeHgla), ggdraw(ptreeMR))
pbis <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 6))
pbis
ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "_with_tree_fixed_params_noBM_effect_size_", effect_size_plot, ".pdf")),
       height = 1.3*columnwidth,
       width = columnwidth,
       unit = "in")

###############################################################################
## Figure S3: design and methods comparisons 
################################################################################
effect_size_plot <- 4
df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM")
             & use_lengths ==  "with_lengths"
             & tree_type == "three_rep"
             & transformation %in% c("log2")
             & effect_size == effect_size_plot
             & fact_disp == base_fact_disp
             & score %in% plot_scores
)
df <- df[!(grepl("_raw", df$methodMod)), ]
df <- df[!(grepl("notrend", df$methodMod)), ]
df <- df[!(grepl("noeBayes", df$methodMod)), ]
df <- df[!(grepl("noreg", df$methodMod)), ]
df <- df[!is.na(df$methodMod), ]
df <- subset(df, fitprocess != "BM")
df <- subset(df, model_process != "BM")

## plot results
p <- ggplot(df, aes(x = emp, y = score_value, fill = methodBare, color = methodBare)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = TRUE,
               position = "identity", width = 0.2,
               aes(group = interaction(methodMod, score, emp, model_process, Design, tree_type))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, use_lengths, tree_type, model_process, methodMod, Design), colour = methodBare)) +
  scale_fill_manual(values = color_methods, name = "") +
  scale_color_manual(values = color_methods,  name = "") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.5, -0.15),
        legend.justification = c("center", "top"),
        legend.key.size = unit(10, 'pt'),
        legend.direction = "horizontal",
        # legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.5,-0.5), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
p

## add tree
bottom_row <- p + theme(plot.margin = unit(c(-0.3,0.1,0.5,-0.5), "cm"))
top_row <- plot_grid(ggdraw(ptreeHgla), ggdraw(ptreeMR))
pbis <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 6))
pbis
ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "_with_tree_all_sim_process_noBM_design_three_rep_effectsize_", effect_size_plot, ".pdf")),
       height = 1.3*columnwidth,
       width = columnwidth,
       unit = "in")

###############################################################################
## Figure S4: eBayes no eBayes trend no trend
################################################################################
effect_size_plot <- 4
df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM")
             & use_lengths ==  "with_lengths"
             & tree_type == "three_rep"
             & transformation %in% c("log2")
             & effect_size == effect_size_plot
             & fact_disp == base_fact_disp
             & score %in% plot_scores
)
df <- df[(grepl("phyloDE", df$methodMod)), ]
df <- subset(df, fitprocess != "BM")
df <- df[!(grepl("noreg", df$methodMod)), ]
df <- df[!is.na(df$methodMod), ]

df$ebtrend <- "no EB"
df$ebtrend[df$eBayes & df$trend == "with_trend"] <- "with trend"
df$ebtrend[df$eBayes & df$trend == "no_trend"] <- "no trend"
df$ebtrend <- factor(df$ebtrend, levels = c("no EB", "no trend", "with trend"))

df <- subset(df, fitprocess != "BM")
df <- subset(df, model_process != "BM")
df <- subset(df, ebtrend != "no EB")

if (nrow(df) / 3 / 2 / 2 / 2 != Nmax - Nmin + 1) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

dodge <- position_dodge(width = 0)

## plot results
p <- ggplot(df, aes(x = emp, y = score_value, fill = ebtrend, color = ebtrend)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(rows = vars(score), cols = vars(Design), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = TRUE,
               position = "identity", width = 0.2,
               aes(group = interaction(methodMod, score, emp, model_process, Design, tree_type))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, use_lengths, tree_type, model_process, methodMod, Design), colour = ebtrend)) +
  scale_fill_manual(values = color_ebayes[c(1, 3)], name = "") +
  scale_colour_manual(values = color_ebayes[c(1, 3)], name = "") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.5, -0.15),
        legend.justification = c("center", "top"),
        legend.key.size = unit(10, 'pt'),
        legend.direction = "horizontal",
        # legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.5,-0.5), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
p

## add tree
bottom_row <- p + theme(plot.margin = unit(c(-0.3,0.1,0.5,-0.5), "cm"))
top_row <- plot_grid(ggdraw(ptreeHgla), ggdraw(ptreeMR))
pbis <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 6))
pbis
ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "_with_tree_all_sim_process_noBM_design_three_rep_eBayes_effectsize_", effect_size_plot, ".pdf")),
       height = 1.3*columnwidth,
       width = columnwidth,
       unit = "in")

###############################################################################
## Figure S5: design and methods comparisons 
################################################################################
effect_size_plot <- 4
df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM")
             & use_lengths ==  "with_lengths"
             & tree_type == "three_rep"
             & transformation %in% c("log2")
             & effect_size == effect_size_plot
             & fact_disp == base_fact_disp
             & score %in% plot_scores
)
df <- df[!(grepl("_raw", df$methodMod)), ]
df <- df[!(grepl("notrend", df$methodMod)), ]
df <- df[!(grepl("noeBayes", df$methodMod)), ]
df <- df[!(grepl("noreg", df$methodMod)), ]
df <- df[(grepl("phyloDE", df$methodMod)), ]
df <- df[!is.na(df$methodMod), ]

if (nrow(df) / 4 / 2 / 3 / 2 != Nmax - Nmin + 1) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

df$fitprocess <- factor(df$fitprocess, levels = c("BM", "OU"))
levels(df$fitprocess) <- c("phyloDE BM", "phyloDE OU")
dodge <- position_dodge(width = 0)

ggplot(df, aes(x = emp, y = score_value, fill = fitprocess, color = fitprocess)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(rows = vars(score), cols = vars(Design, model_process), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = TRUE,
               position = "identity", width = 0.2,
               aes(group = interaction(methodMod, score, emp, model_process, Design, tree_type))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(effect_size, use_lengths, tree_type, model_process, methodMod, Design), colour = fitprocess)) +
  scale_linetype(name = "") +
  scale_fill_manual(values = color_methods[c(1, 4)], name = "") +
  scale_color_manual(values = color_methods[c(1, 4)], name = "") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.5, -0.15),
        legend.justification = c("center", "top"),
        legend.key.size = unit(10, 'pt'),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        # legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,0.2,-0.5), "cm"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
# colorblindr::cvd_grid()
ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "_all_sim_process_design_effectsize_", effect_size_plot, ".pdf")),
       height = 1.2 * columnwidth,
       width = 2 * columnwidth,
       unit = "in")

###############################################################################
## Figure S1: fixed parameters effect size varying
################################################################################
df <- subset(result.table_plot,
             lengthNormalization %in% c("TPM")
             & use_lengths ==  "with_lengths"
             & transformation %in% c("log2")
             & emp == "fixed"
             & score %in% plot_scores
)
df <- df[!(grepl("_raw", df$methodMod)), ]
df <- df[!(grepl("notrend", df$methodMod)), ]
df <- df[!(grepl("noeBayes", df$methodMod)), ]
df <- df[!(grepl("noreg", df$methodMod)), ]
df <- df[!is.na(df$methodMod), ]
df <- subset(df, fitprocess != "BM")
df <- subset(df, model_process != "BM")

if (nrow(df) / 4 / 2 / 4 / 3 != Nmax - Nmin + 1) message("wrong dimension")
if (anyNA(df$score_value)) message("Some NAs in the data")

dodge <- position_dodge(width = 0)

df$tree_type <- as.factor(df$tree_type)
levels(df$tree_type) <- c("original tree", "three replicates")
ggplot(df, aes(x = effect_size, y = score_value, fill = methodBare, color = methodBare)) +
  geom_hline(aes(yintercept = threshold), linetype = 2) +
  facet_grid(rows = vars(score), cols = vars(Design, tree_type), scales = "free", switch = "y", labeller = labeller(Design = label_both)) +
  geom_boxplot(alpha = 0.3, outlier.size = 0.5, show.legend = TRUE,
               position = "identity", width = 0.2,
               aes(group = interaction(methodMod, score, emp, model_process, Design, tree_type, effect_size))) +
  stat_summary(fun = median, geom = 'line', aes(group = interaction(use_lengths, tree_type, model_process, methodMod, Design), colour = methodBare)) +
  scale_fill_manual(values = color_methods, name = "") +
  scale_color_manual(values = color_methods,  name = "") +
  xlab("Effect Size") +
  ylab("") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.5, -0.15),
        legend.justification = c("center", "top"),
        legend.key.size = unit(10, 'pt'),
        legend.direction = "horizontal",
        plot.margin = unit(c(0.1,0.1,0.7,-0.5), "cm"),
  )

ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "_fixed_params_varying_effect_size.pdf")),
       height = columnwidth,
       width = 2 * columnwidth,
       unit = "in")

###############################################################################
## Figure S2: distribution of simulation parameters
################################################################################

df_params <- rbind(
  data.frame(parameter = "alpha",
             var = "low var",
             value = alphas_OU_sim_low),
  data.frame(parameter = "alpha",
             var = "high var",
             value = alphas_OU_sim_high),
  data.frame(parameter = "lambda",
             var = "low var",
             value = lambdas_OU_sim_low),
  data.frame(parameter = "lambda",
             var = "high var",
             value = lambdas_OU_sim_high)
)
df_params$var <- factor(df_params$var, levels = c("low var", "high var"))

df_median <- rbind(
  data.frame(parameter = "alpha",
             value = alpha_median),
  data.frame(parameter = "lambda",
             value = lambda_median)
)

ggplot(df_params, aes(x=value, color = var, fill = var)) +
  facet_wrap(vars(parameter), scales = "free", strip.position = "bottom") +
  geom_vline(data = df_median, aes(xintercept = value), linetype = "dashed") + 
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha = 0.5) + 
  scale_fill_manual(values = color_methods[c(1, 4)], name = "") +
  scale_color_manual(values = color_methods[c(1, 4)], name = "") +
  xlab("") +
  ylab("density") +
  theme_bw() +
  theme(text = element_text(size = 9),
        title = element_text(size = 7),
        strip.text.y = element_text(vjust = 0),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.53, 0.99),
        legend.justification = c("left", "top"),
        legend.key.size = unit(10, 'pt'),
        # legend.position = "bottom",
        plot.margin = unit(c(0.1,0.1,-0.3,0.1), "cm"),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

# colorblindr::cvd_grid()
ggsave(file.path(results_directory, paste0(datestamp_day_anaysis, "distribution_parameters.pdf")),
       height = 0.6*columnwidth,
       width = 1.5 * columnwidth,
       unit = "in")
