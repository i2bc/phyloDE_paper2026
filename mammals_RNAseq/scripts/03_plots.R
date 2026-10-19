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
# Saved in `mammals_RNAseq/results/YYYY-MM-DD_YYYY-MM-DD_simulations_results_mammals_RNAseq`
################################################################################

## OUTPUT ######################################################################
# Analysis plots
################################################################################

library(ggplot2)
library(here)
library(RColorBrewer)
library(ggtree)
library(patchwork)
library(viridisLite)
library(ape)

###############################################################################
# Standardize color palette for all figures in the manuscript
###############################################################################

color_methods    <- mako(4, alpha = 1, begin = 0.1, end = 0.8, direction = 1)[1:4]
color_methods[4] <- inferno(4, alpha = 1, begin = 0, end = 0.8, direction = 1)[4]
color_palete     <- scale_fill_manual(values = color_methods[c(1,3,2,4)], guide = "none")
fill_palete      <- scale_color_manual(values = color_methods[c(1,3,2,4)], guide = "none")
columnwidth      <- 3

###############################################################################
# Load results and subset to OU models for phyloDE and phylolm
###############################################################################

here_dir        <- here("mammals_RNAseq", "results")
datestamp_simu  <- "2026-02-09"
simus_directory <- file.path(here_dir, paste0(datestamp_simu, "_simulations_mammals"))
designs         <- readRDS(file.path(simus_directory, "designs.rds"))
datestamp_res   <- "2026-02-09"
res_directory   <- file.path(here_dir, paste0(datestamp_simu, "_", datestamp_res, "_simulations_results_mammals"))
res             <- readRDS(file.path(res_directory, "res.rds"))
res$row_id      <- 1:nrow(res)
res[,c("MCC", "TP", "FP", "TN", "FN")] <- NULL
res             <- res[res$method %in% c("phylolmOU", "eve", "limma", "phyloDEOU"),]
res$method      <- as.character(res$method)
res$method[res$method=="phylolmOU"] <- "phylolm"
res$method[res$method=="phyloDEOU"] <- "phyloDE"

###############################################################################
# Generate Plot
###############################################################################

for (th in c(0.05, 0.1)) {
  res_long <- reshape(res[res$th==th,],
                      varying   = c("FDR", "TPR"),
                      v.names   = "value",
                      timevar   = "metric",
                      times     = c("FDR", "TPR"),
                      idvar     = "row_id",
                      direction = "long")
  
  res_long$metric <- factor(res_long$metric, c("FDR", "TPR"))
  
  res_long$method <- factor(res_long$method, c("limma", "phylolm", "eve", "phyloDE"))
  
  tree_plots <- lapply(c("primates","rodents"), function(k) {
    
    tree           <- designs$treeTr
    species_cond2  <- names(designs$allConds[[k]])[designs$allConds[[k]]==2]
    tree           <- groupClade(tree, getMRCA(tree, species_cond2))
    musind <- grepl("mus.*", tree$tip.label)
    tree$tip.label <- sub("mus", "m", tree$tip.label)
    tree$tip.label[!musind] <- substring(tree$tip.label[!musind], 1, 3)
    tree$tip.label[musind] <- tolower(substring(tree$tip.label[musind], 1, 4))
    # tree$tip.label[duplicated(tree$tip.label)] <- ""
    
    ggtree(tree, aes(color=factor(group)), size=0.5) +
      geom_tiplab(size=2.5, hjust=0, offset = 0.01) +
      scale_color_manual(values = c("grey65", "black")) +
      theme_tree() +
      guides(color = "none") +
      xlim(c(0, 1.2))
    # ggtitle(paste("Effective sample size: ", round(designs$nEff[[k]], 2))) +
    # theme(plot.title = element_text(hjust = 0.5, size=30))
    
  })
  
  main_plot <- ggplot(res_long, aes(x=method, y=value, color=method, fill=method)) +
    facet_grid(metric~condType, switch = "y") +
    geom_boxplot(alpha = 0.3, outlier.size = 0.5) +
    geom_hline(data = data.frame(metric=factor("FDR", c("FDR", "TPR")), y=th), 
               aes(yintercept = y), 
               linetype = 2) +
    guides(color = "none") +
    color_palete + 
    fill_palete +
    ylab("") +
    theme_bw() +
    theme(text             = element_text(size=9), 
          strip.placement  = "outside",
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x     = element_blank(),
          strip.text.y     = element_text(vjust = 0),
          axis.title.x     = element_text(margin = margin(t = 10)),
          axis.title.y     = element_text(margin = margin(r = 0)),
    )
  
  
  pdf(width=2*columnwidth, height=2.5*columnwidth, file.path(res_directory, paste0("plot_stats_",th,".pdf")))
  print(wrap_plots(tree_plots, ncol=2) / main_plot +
          plot_layout(heights = c(5, 5)))
  dev.off()
}
