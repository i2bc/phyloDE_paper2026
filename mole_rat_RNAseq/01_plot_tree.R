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
# Program to plot the mole rat tree, 
# transformed for several values of rho and lambda.
################################################################################

## INPUT #######################################################################
# Formatted data as a compcodeR object
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Figure 3
# Saved in `mole_rat_RNAseq/data/data_formatted`
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
plot(tree_3rep)


################################################################################
## Utility functions
################################################################################
trans_alpha <- function(tree, alpha, lambda) {
  tree_model <- phylolm::transf.branch.lengths(tree, "OUfixedRoot", parameters = list(alpha = alpha))$tree
  tree_model <- phylolm::transf.branch.lengths(tree_model, "lambda", parameters = list(lambda = lambda))$tree
  return(tree_model)
}

get_rho <- function(alpha, t_tree = 1, tol = .Machine$double.eps) {
  zeros <- alpha <= tol
  nas <- is.na(alpha)
  res <- rep(1, length(alpha))
  res[!zeros & !nas] <- 1 - -expm1(-2 * alpha[!zeros & !nas] * t_tree) / (2 * alpha[!zeros & !nas] * t_tree)
  res[nas] <- NA
  return(res)
}

get_rho_inv <- function(r, t_tree = 1, alpha_bounds = c(0.0001, 300)) {
  if (r == 0) return(0)
  uniroot((function(x) get_rho(x, t_tree) - r), interval = alpha_bounds, tol = .Machine$double.eps^0.5)$root
}

################################################################################
## Plot
################################################################################
columnwidth <- 3 # in inches # 1 point / 72 = 1 inch
results_directory <- here(data_file, "data_formatted")

lambda_vals <- c(0.1, 0.25, 0.5, 0.75, 0.9) # c(0.75, 0.5, 0.25)
rho_vals <- c(0.13, 0.45, 0.66, 0.96) # c(0.25, 0.5, 0.75)
nrow <- length(rho_vals) + 2
ncol <- length(lambda_vals) + 2

pdf(file.path(results_directory, "tree_transformations.pdf"),
    height = 1.3*20/26*columnwidth,
    width = 1.3*columnwidth)
layout(mat = matrix(1:(ncol*nrow), nrow = nrow, ncol = ncol),
       heights = c(3, 1, rep(5, nrow - 2)),
       widths = c(3, rep(5, ncol - 2), 1.5)) 
# layout.show(16)
par(mar = c(0.1, 0.1, 0.1, 0.1) + 0.1)
# B
plot.new()
plot.new()
for (rho in rho_vals) {
  plot.new()
  if (rho == 0.66) {
    plot.window(c(0, 1), c(0, 1))
    text(0.2, 0.5, "B", font = 2, cex = 1.3)
    arrows(0.5, 0.5, 0.9, 0.5, lwd = 2, length = 0.05)
  }
}
# trees
for (lambda in lambda_vals) {
  plot.new()
  if (lambda == 0.75) {
    plot.window(c(0, 1), c(0, 1))
    text(0.5, 0.8, "A", font = 2, cex = 1.3)
    arrows(0.5, 0.5, 0.5, 0.1, lwd = 2, length = 0.05, code = 2)
  }
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  text(0.5, 0.5, bquote(lambda==.(lambda)))
  for (rho in rho_vals) {
    treetrans <- trans_alpha(tree_3rep, alpha = get_rho_inv(rho), lambda = lambda)
    plot(treetrans, show.tip.label = FALSE)
    if (lambda == 0.1 && rho == 0.66) {
      par(xpd = NA)
      rect(-0.1, 0, 6.1, 13, lwd = 2.5)
    }
    if (lambda == 0.75 && rho == 0.96) {
      par(xpd = NA)
      rect(-0.1, 0, 1.1, 55, lwd = 2.5)
    }
  }
}
# labels
plot.new()
plot.new()
for (rho in rho_vals) {
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  text(0.4, 0.5, bquote(rho==.(rho)), srt = -90)
}
dev.off()
