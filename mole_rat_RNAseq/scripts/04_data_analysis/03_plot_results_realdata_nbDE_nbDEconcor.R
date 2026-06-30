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
# Program to perform CRE enrichment analyses
################################################################################

## INPUT #######################################################################
# Results from script `02_lfc_cre_vs_lfc_rnaseq.R`
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_data_analysis`
################################################################################

## OUTPUT ######################################################################
# Figure 5 and Tables S1-6
################################################################################

################################################################################
################################################################################

library(ggplot2)
library(dplyr)
library(readr)
library(here)

PATH <- here("mole_rat_RNAseq/results/2026-03-11_data_analysis/")
df <- read_tsv(file.path(PATH, "fischer_lfc_cre_vs_lfc_rnaseq_by_tissue_branch.tsv"))

df_plot <- df %>%
  filter(
    method %in% c("EVEemp", "phyloDEOU"),
    filter == "adjpvalue < 0.05",
    branch %in% c("MR", "hgla", "fdam"),
    tissue %in% c("Heart", "Liver")
  ) %>%
  mutate(
    Design = recode(branch, "hgla" = "Hgla", "fdam" = "Fdam"),
    Design = factor(Design, levels = c("MR", "Hgla", "Fdam")),
    method = recode(method, "EVEemp" = "eve", "phyloDEOU" = "phyloDE"),
    method = factor(method, levels = c("eve", "phyloDE")),
    tissue = factor(tissue, levels = c("Heart", "Liver")),
    rate_label = paste0(round(rate * 100), "%")
  )

gg <- ggplot(df_plot, aes(x = method)) +
  geom_col(
    aes(y = n, fill = "Total n"),
    width = 0.7,
    alpha = 0.45
  ) +
  geom_col(
    aes(y = concordant, fill = "Concordant"),
    width = 0.7
  ) +
  geom_text(
    aes(y = n, label = rate_label),
    vjust = -0.4,
    size = 2.8
  ) +
  facet_grid(tissue ~ Design, labeller = labeller(Design = label_both)) +
  labs(
    x = element_blank(),
    y = "Number of DE genes",
    fill = "",
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 200)) +
  theme(
    text = element_text(size = 9),
    title = element_text(size = 7),
    # strip.text.y = element_text(vjust = 0),
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(10, 'pt'),
    legend.margin = margin(t = -4),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
  ) +
  scale_fill_manual(
    values = c(
      "Total n" = "#4F8EAD",
      "Concordant" = "#1f78b4"
    ),
    labels = c(
      "Total n" = "Total DE",
      "Concordant" = "Concordant"
    )
  )

gg
columnwidth <- 3

ggsave(file.path(PATH, "realdata_heart_liver_concordant.pdf"),
       height = columnwidth,
       width = columnwidth,
       unit = "in")

###############################################################################
## Plot tree
################################################################################
library(ape)
library(cowplot)

# load tree and drop replicates
tree <- read.tree(file = here("mole_rat_RNAseq/data/data_formatted/species_tree_mya.nwk"))
tree <- drop.tip(tree, tree$tip.label[duplicated(sub("_.*$", "", tree$tip.label))])
tree$edge.length <- tree$edge.length / max(vcv(tree))
tree$tip.label <- sub("_.*$", "", tree$tip.label)

col_orig <- c("grey60", "black")

ptreeMR <- function() {
  par(mar = c(0, 2, 0, 0))
  col_plot <-  rep(col_orig[1], 6)
  col_plot[c(2, 3, 4)] <- col_orig[2]
  col_tip <-  rep(col_orig[1], 4)
  col_tip[c(1, 2)] <- col_orig[2]
  plot(tree, edge.color = col_plot, edge.width = 3,
       y.lim = c(-0.1, 4.5), x.lim = c(0, 2),
       show.tip.label = TRUE, cex = 0.6, label.offset = 0.1,
       tip.color = col_tip)
}
# ggdraw(ptreeMR)

ptreeHgla <- function() {
  par(mar = c(0, 1, 0, 1))
  col_plot <-  rep(col_orig[1], 6)
  col_plot[3] <- col_orig[2]
  col_tip <-  rep(col_orig[1], 4)
  col_tip[1] <- col_orig[2]
  plot(tree, edge.color = col_plot, edge.width = 3,
       y.lim = c(-0.1, 4.5), x.lim = c(0, 2),
       show.tip.label = TRUE, cex = 0.6, label.offset = 0.1,
       tip.color = col_tip)
  #
}
# ggdraw(ptreeHgla)

ptreeFdam <- function() {
  par(mar = c(0, 0, 0, 2))
  col_plot <-  rep(col_orig[1], 6)
  col_plot[4] <- col_orig[2]
  col_tip <-  rep(col_orig[1], 4)
  col_tip[2] <- col_orig[2]
  plot(tree, edge.color = col_plot, edge.width = 3,
       y.lim = c(-0.1, 4.5), x.lim = c(0, 2),
       show.tip.label = TRUE, cex = 0.6, label.offset = 0.1,
       tip.color = col_tip)
  #
}
# ggdraw(ptreeFdam)

## with tree
top_row <- plot_grid(ggdraw(ptreeMR), ggdraw(ptreeHgla), ggdraw(ptreeFdam), nrow = 1)
pbis1 <- gg + theme(plot.margin = unit(c(-0.3,0.1,0.1,0.1), "cm"))

pbis1 <- plot_grid(top_row, pbis1, ncol = 1, rel_heights = c(1, 6))
pbis1

ggsave(file.path(PATH, "realdata_heart_liver_concordant_with_tree.pdf"),
       height = 1.1*columnwidth,
       width = columnwidth,
       unit = "in")

###############################################################################
## Gene list
################################################################################
library(kableExtra)

rnaseq_daunesse <- read.table(here("mole_rat_RNAseq/data/Daunesse2025/combined_deseq2Xeve_results.tab"), header=TRUE)

de_genes <- read.table(file.path(PATH, "mole_rat_RNAseq_data.tsv"))
de_genes <- subset(de_genes, method == "phyloDEOU" & adjpvalue < 0.05)
tables <- NULL

for (tissue_sub in c("Liver", "Heart")) {
  for (branch_sub in c("MR", "hgla", "fdam")) {
    de_genes_sub <- subset(de_genes, tissue == tissue_sub & branch == branch_sub)
    de_genes_sub$gene.name <- rnaseq_daunesse$GeneName[match(de_genes_sub$gene.id, rnaseq_daunesse$gene.id)]
    de_genes_sub <- de_genes_sub[, c("gene.id", "gene.name", "baseMean", "logFC", "adjpvalue")]
    de_genes_sub <- de_genes_sub[order(de_genes_sub$adjpvalue, -abs(de_genes_sub$logFC)), ]
    de_genes_sub$adjpvalue <- de_genes_sub$adjpvalue * 100
    rownames(de_genes_sub) <- NULL
    colnames(de_genes_sub) <- c("gene id", "gene name", "base log mean", "log fold change", "adjusted p-value (\\%)")
    
    tables <- de_genes_sub %>%
      kable("latex", digits = 2, escape = F, align = "r", booktabs = T, longtable = T, linesep = "",
            caption = paste0(nrow(de_genes_sub), " DE genes found by \\code{phyloDE} for the ", tissue_sub, " data with the ``", branch_sub, "'' design.")) %>%
      # column_spec(c(5, 9), bold = TRUE) %>%
      # add_header_above(c(" ", "ML" = 4, "REML" = 4)) %>%
      # column_spec(5, border_right = T) %>%
      # row_spec(2, hline_after = T) %>%
      # row_spec(5, hline_after = T) %>%
      # group_rows("", 3, 5, latex_gap_space = "0em") %>%
      kable_styling(latex_options = c("hold_position", "repeat_header"), position = "center") %>%
      kable_styling(font_size = 8)
    writeLines(tables, file.path(PATH, paste0("mole_rat_RNAseq_", tissue_sub, "_", branch_sub, ".tex")))
    write_tsv(de_genes_sub, file.path(PATH, paste0("mole_rat_RNAseq_", tissue_sub, "_", branch_sub, ".tsv")))
  }
}
