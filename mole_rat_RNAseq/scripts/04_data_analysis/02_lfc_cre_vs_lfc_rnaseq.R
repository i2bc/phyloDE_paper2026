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
# Results of `phyloDE` and `evemodel`
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_data_analysis`
################################################################################

## OUTPUT ######################################################################
# Analysis of CRE enrichment of lists of genes
################################################################################

library(readxl)
library(data.table)
library(ggplot2)
library(here)

################################################################################
## Parameters
################################################################################
datestamp_day <- "2026-03-11" #format(Sys.time(), "%Y-%m-%d")
data_file <- "mole_rat_RNAseq/data"
results_dir <- here("mole_rat_RNAseq", "results", paste0(datestamp_day, "_data_analysis"))

sp      <- c("Hgla", "Fdam", "Cpor", "Mmus")
designs <- list(hgla = c(2,1,1,1),
                MR   = c(2,2,1,1),
                fdam = c(1,2,1,1))

# Each filter defines a (method, column, threshold) combination matching one scatter plot
filters <- list(
  list(method = "EVEDaunesse",    col = "pvalue",    thresh = 0.05),
  #list(method = "EVEDaunesse",    col = "adjpvalue", thresh = 0.05),
  list(method = "EVEemp",         col = "pvalue",    thresh = 0.05),
  list(method = "EVEemp",         col = "adjpvalue", thresh = 0.05),
  list(method = "phyloDEOU",      col = "adjpvalue", thresh = 0.05),
  list(method = "phyloDEBM",      col = "adjpvalue", thresh = 0.05),
  list(method = "phyloDEOUbatch", col = "adjpvalue", thresh = 0.05),
  list(method = "phyloDEBMbatch", col = "adjpvalue", thresh = 0.05)
)

################################################################################
## Load data
################################################################################
rnaseq <- as.data.table(rbind(read.table(here(results_dir, "mole_rat_RNAseq_data.tsv")),
                              read.table(here(results_dir, "mole_rat_RNAseq_data_batch.tsv"))))

rnaseq_daunesse <- as.data.table(read.table(here("mole_rat_RNAseq/data/Daunesse2025/combined_deseq2Xeve_results.tab"), 
                                            header=TRUE))  

#rnaseq_daunesse$baseMean  <- NA_real_
rnaseq_daunesse[branch == "Fdam", branch := "fdam"]                                                                           
rnaseq_daunesse[branch == "Hgla", branch := "hgla"]
rnaseq_daunesse[, adjpvalue := p.adjust(empirical.pval, method = "BH"), by = .(tissue, branch)] # none pass the threshold if pval corrected
rnaseq_daunesse$method <- "EVEDaunesse"
rnaseq_daunesse <- rnaseq_daunesse[, .(gene.id,  tissue, branch, pvalue = empirical.pval, adjpvalue, logFC = shift, baseMean, method)]
rnaseq <- rbind(rnaseq, rnaseq_daunesse)

cre    <- as.data.table(
  read_excel(here("mole_rat_RNAseq/data/Daunesse2025/tables 2/S5_regulatory_score.xlsx"), sheet = 3)
)

################################################################################
## Compute CRE logFC per gene x tissue x branch (data.table, vectorised)
##
## Strategy: expand cre by joining a species-to-group mapping table, then
## compute group means in one grouped aggregation
################################################################################
sp_groups <- rbindlist(lapply(names(designs), function(b) {
  data.table(branch = b, species = sp, group = designs[[b]])
}))

# Join expands each (gene, tissue, species) row into one row per branch
cre_tr <- cre[sp_groups, on = "species", allow.cartesian = TRUE][
  , .(logFC.cre = log2(mean(reg_score_quant_norm[group == 2] + 0.01)) -
                  log2(mean(reg_score_quant_norm[group == 1] + 0.01))),
  by = .(gene.id = gene_id, tissue, branch)
]

################################################################################
## Merge CRE logFC with RNAseq differential expression results
################################################################################
cre_X_rnaseq <- merge(cre_tr, rnaseq, by = c("gene.id", "tissue", "branch"))
setorder(cre_X_rnaseq, gene.id, tissue, branch)

################################################################################
## Fisher test on sign concordance (sign of logFC.rna vs sign of logFC.cre)
##
## 2x2 table:  sign(logFC.rna) x sign(logFC.cre)
################################################################################
fisher_sign <- rbindlist(lapply(filters, function(f) {
  sub <- cre_X_rnaseq[method == f$method & get(f$col) < f$thresh]
  tab <- table(sign_rna = sign(sub$logFC), sign_cre = sign(sub$logFC.cre))
  tab <- tab[rownames(tab) %in% c("-1","1"), colnames(tab) %in% c("-1","1")]
  ft  <- fisher.test(tab)
  n   <- nrow(sub)
  con <- sum(sign(sub$logFC) == sign(sub$logFC.cre), na.rm = TRUE)
  data.table(method     = f$method,
             filter     = paste0(f$col, " < ", f$thresh),
             n          = n,
             concordant = con,
             rate       = round(con / n, 2),
             OR         = round(ft$estimate, 2),
             pvalue     = ft$p.value)
}))

print(fisher_sign, row.names = FALSE)

write.table(fisher_sign,
            here(results_dir, "fischer_lfc_cre_vs_lfc_rnaseq.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

################################################################################
## Summary figure: 4-panel scatter (logFC.rna vs logFC.cre) with annotations
################################################################################
# Build a combined dataset with a panel label column
panel_label <- function(f) paste0(f$method, "\n", f$col, " < ", f$thresh)

plot_data <- rbindlist(lapply(filters, function(f) {
  sub        <- cre_X_rnaseq[method == f$method & get(f$col) < f$thresh]
  sub$label  <- panel_label(f)
  sub
}))
plot_data$label <- factor(plot_data$label,
                          levels = sapply(filters, panel_label))

p <- ggplot(plot_data, aes(x = logFC, y = logFC.cre)) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(x = expression(log[2]*FC~(RNA-seq)),
       y = expression(log[2]*FC~(CRE~score))) +
  theme_bw(base_size = 20)
p

ggsave(here(results_dir, "lfc_cre_vs_lfc_rnaseq.pdf"),
       plot = p, width = 14, height = 14)

################################################################################
## Same fisher test and scatter plot, split by (tissue, branch)
################################################################################
tissue_branch <- unique(cre_X_rnaseq[, .(tissue, branch)])

fisher_sign_by_tb <- rbindlist(lapply(filters, function(f) {
  rbindlist(lapply(seq_len(nrow(tissue_branch)), function(i) {
    ti <- tissue_branch$tissue[i]
    br <- tissue_branch$branch[i]
    sub <- cre_X_rnaseq[method == f$method & get(f$col) < f$thresh &
                          tissue == ti & branch == br]
    if (nrow(sub) < 2) return(NULL)
    tab <- table(sign_rna = sign(sub$logFC), sign_cre = sign(sub$logFC.cre))
    tab <- tab[rownames(tab) %in% c("-1","1"), colnames(tab) %in% c("-1","1")]
    if (any(dim(tab) < 2)) return(NULL)
    ft  <- fisher.test(tab)
    n   <- nrow(sub)
    con <- sum(sign(sub$logFC) == sign(sub$logFC.cre), na.rm = TRUE)
    data.table(tissue     = ti,
               branch     = br,
               method     = f$method,
               filter     = paste0(f$col, " < ", f$thresh),
               n          = n,
               concordant = con,
               rate       = round(con / n, 2),
               OR         = round(ft$estimate, 2),
               pvalue     = ft$p.value)
  }))
}))

print(fisher_sign_by_tb, row.names = FALSE)

write.table(fisher_sign_by_tb,
            here(results_dir, "fischer_lfc_cre_vs_lfc_rnaseq_by_tissue_branch.tsv"),
            sep = "\t", quote = TRUE, row.names = FALSE)

plot_data_tb <- rbindlist(lapply(filters, function(f) {
  sub       <- cre_X_rnaseq[method == f$method & get(f$col) < f$thresh]
  sub$label <- panel_label(f)
  sub
}))
plot_data_tb$label <- factor(plot_data_tb$label,
                             levels = sapply(filters, panel_label))

p_tb <- ggplot(plot_data_tb, aes(x = logFC, y = logFC.cre)) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  facet_grid(tissue + branch ~ label, scales = "free") +
  labs(x = expression(log[2]*FC~(RNA-seq)),
       y = expression(log[2]*FC~(CRE~score))) +
  theme_bw(base_size = 12)
p_tb

n_panels_x <- length(filters)
n_panels_y <- nrow(tissue_branch)

ggsave(here(results_dir, "lfc_cre_vs_lfc_rnaseq_by_tissue_branch.pdf"),
       plot = p_tb, width = 4 * n_panels_x, height = 4 * n_panels_y)

