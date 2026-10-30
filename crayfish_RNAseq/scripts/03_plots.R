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
# Program to plot the formatted results
################################################################################

## INPUT #######################################################################
# Formatted result table
# Produced by crayfish_RNAseq/scripts/02_analysis.R
# Saved in `crayfish_RNAseq/results/YYYY-MM-DD_YYYY-MM-DD_simulations_results_crayfish_RNAseq`
# The full result table from our own run is provided for convenience
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


here_dir        <- here("crayfish_RNAseq", "results")
datestamp_simu  <- "2026-01-19"
simus_directory <- file.path(here_dir, paste0(datestamp_simu, "_simulations_crayfish"))
designs         <- readRDS(file.path(simus_directory, "designs.rds"))
datestamp_res   <- "2026-01-19"
res_directory   <- file.path(here_dir, paste0(datestamp_simu, "_", datestamp_res, "_simulations_results_crayfish"))
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
  
  tree_plots <- lapply(c(2,4,6), function(k) {
    
    tree           <- designs[[1]][[1]][[k]]$treeTr
    species_cond2  <- substring(designs[[1]][[1]][[k]]$species_cond2, 1, 2)
    tree$tip.label <- substring(tree$tip.label, 1, 2)
    tree           <- groupClade(tree, getMRCA(tree, species_cond2))
    tree$tip.label[duplicated(tree$tip.label)] <- ""
    
    ggtree(tree, aes(color=factor(group)), size=0.5) +
      geom_tiplab(size=2.5, hjust=-0.3) +
      scale_color_manual(values = c("grey65", "black")) +
      theme_tree() +
      guides(color = "none") +
      xlim(c(0, 1.5))
    # ggtitle(paste("Effective sample size: ", round(designs[[1]][[1]][[k]]$neff, 2))) +
    # theme(plot.title = element_text(hjust = 0.5, size=30))
  })
  
  main_plot <- ggplot(res_long, aes(x=method, y=value, color=method, fill=method)) +
    facet_grid(metric~design, switch = "y") +
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
  
  pdf(width=2*columnwidth, height=2*columnwidth, file.path(res_directory, paste0("plot_stats_",th,".pdf")))
  print(wrap_plots(tree_plots, ncol=3) / main_plot +
          plot_layout(heights = c(2, 3)))
  dev.off()
  
}
