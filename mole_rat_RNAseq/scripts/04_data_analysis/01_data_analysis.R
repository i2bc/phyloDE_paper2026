# Copyright (C) {2025} {PB, MG}
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
# Program to analyse mole rat data
################################################################################

## INPUT #######################################################################
# Formatted data
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################

## OUTPUT ######################################################################
# Results of `phyloDE` and `evemodel`
# Saved in `mole_rat_RNAseq/results/YYY-MM-DD_data_analysis`
################################################################################

library(evemodel)
library(ape)
library(phytools)
library(ggplot2)
library(here)
library(phyloDE)

################################################################################
## Parameters
################################################################################
datestamp_day <- "2026-03-11" #format(Sys.time(), "%Y-%m-%d")
data_file <- "mole_rat_RNAseq/data/data_formatted"
results_directory <- here("mole_rat_RNAseq", "results", paste0(datestamp_day, "_data_analysis"))
dir.create(results_directory, showWarnings = FALSE)
data_name <- "mole_rat_RNAseq_data"

all_results <- NULL
all_results_batch <- NULL

for (tissue in c("Liver", "Heart")) {
  for (branch in c("hgla", "MR", "fdam")) {
    ################################################################################
    ## Load data
    ################################################################################
    # tissue <- "Liver" # one of "Liver" or "Heart"
    # branch <- "hgla"  # one of "MR", "hgla", "fdam" or "cpor"
    
    data_norm <- read.table(here(data_file, paste0(tissue, "_salmon_matrix_normalised_counts.tab")))
    data_norm <- as.matrix(data_norm)
    tree <- read.tree(file = here("mole_rat_RNAseq/data/data_formatted/species_tree_mya.nwk"))
    
    # problem with many zero entries in the reg matrix ??
    sum(apply(data_norm, 1, function(xx) any(xx == 0)))
    mean(apply(data_norm, 1, function(xx) any(xx == 0))) # 29% of genes have at least one 0
    # heatmap(data_norm)
    
    # Filtering genes with zero variance
    # See https://gitlab.com/maelle.daunesse/daunesse_et_al_2025/-/blob/main/scripts/02_transcriptomic_analyses/03_shift_analysis/01_EVE_analysis.R?ref_type=heads
    zeroVarRows <- apply(data_norm, 1, function(x) var(x) < 1e-6)
    sum(zeroVarRows)
    data_norm <- data_norm[!zeroVarRows, ]
    
    ################################################################################
    ## Design
    ################################################################################
    
    # Function to get parent branches of the set of tips.
    # tips must be named 
    get_parent_branches <- function(phy, tips) {
      # if only one tip, just get the parent branch
      if (length(tips) == 1) {
        return((phy$edge[, 2] == which(phy$tip.label == tips)) + 0)
      }
      # if several tips, get the MRCA, then descendants
      par <- getMRCA(phy, tips)
      par <- c(par, getDescendants(phy, par))
      return(phy$edge[, 2] %in% par + 0)
    }
    
    ## design for phyloDE and EVE
    if (branch == "MR"){
      condition <- grepl("Hgla|Fdam", colnames(data_norm)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Hgla", "Fdam"))
    } else if (branch == "hgla"){
      condition <- grepl("Hgla", colnames(data_norm)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Hgla"))
    }  else if (branch == "fdam"){
      condition <- grepl("Fdam", colnames(data_norm)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Fdam"))
    }  else if (branch == "cpor"){
      condition <- grepl("Cpor", colnames(data_norm)) + 0
      isTheta2edge <- get_parent_branches(tree, c("Cpor"))
    }
    
    design_formula <- as.formula(~ condition)
    design_data <- data.frame(condition = condition)
    rownames(design_data) <- colnames(data_norm)
    design_data$condition <- factor(design_data$condition)
    design <- model.matrix(design_formula, design_data)
    
    ################################################################################
    ## phyloDE - BM
    ################################################################################

    fitlimmaBM <- phyloDE::phylolmFit(data_norm,
                                         design = design,
                                         phy = tree,
                                         model = 'BM',
                                         measurement_error = TRUE,
                                         use_consensus = TRUE,
                                         ddf_method = "Samples",
                                         medianOU = FALSE,
                                         ncores = 6)

    fitlimmaBM <- limma::eBayes(fitlimmaBM, trend = TRUE)
    limmaBMPvalues <- fitlimmaBM$p.value[, ncol(fitlimmaBM$p.value)]

    resultlimmaBM <- data.frame(gene.id = rownames(data_norm),
                                tissue = tissue,
                                branch = branch,
                                pvalue = limmaBMPvalues,
                                adjpvalue = p.adjust(limmaBMPvalues, method = 'BH'),
                                logFC = fitlimmaBM$coefficients[, ncol(fitlimmaBM$coefficients)],
                                baseMean = fitlimmaBM$coefficients[, 1],
                                method = "phyloDEBM")

    all_results <- rbind(all_results, resultlimmaBM)
    write.table(all_results,
                here(results_directory, paste0(data_name, ".tsv")),
                sep='\t', quote = FALSE)

    ################################################################################
    ## phyloDE - OU
    ################################################################################

    fitlimmaOU <- phyloDE::phylolmFit(data_norm,
                                         design = design,
                                         phy = tree,
                                         model = 'OUfixedRoot',
                                         measurement_error = TRUE,
                                         use_consensus = TRUE,
                                         ddf_method = "Samples",
                                         medianOU = FALSE,
                                         ncores = 6)

    fitlimmaOU <- limma::eBayes(fitlimmaOU, trend = TRUE)
    limmaOUPvalues <- fitlimmaOU$p.value[, ncol(fitlimmaOU$p.value)]

    resultlimmaOU <- data.frame(gene.id = rownames(data_norm),
                                tissue = tissue,
                                branch = branch,
                                pvalue = limmaOUPvalues,
                                adjpvalue = p.adjust(limmaOUPvalues, method = 'BH'),
                                logFC = fitlimmaOU$coefficients[, ncol(fitlimmaOU$coefficients)],
                                baseMean = fitlimmaOU$coefficients[, 1],
                                method = "phyloDEOU")

    all_results <- rbind(all_results, resultlimmaOU)
    write.table(all_results,
                here(results_directory, paste0(data_name, ".tsv")),
                sep='\t', quote = FALSE)

    ################################################################################
    ## EVE - chi2
    ################################################################################
    colSpecies <- sub("_.*$","", colnames(data_norm))

    fitEVE <- twoThetaTest(tree = tree, gene.data = data_norm,
                           isTheta2edge = isTheta2edge, colSpecies = colSpecies,
                           cores = 6)

    ################################################################################
    ## EVE - empirical
    ################################################################################

    # median paramters of the null model estimated from the data
    params <- apply(fitEVE$oneThetaRes$par, 2, function(x) median(x, na.rm = TRUE))

    # run null simulations
    set.seed(1289)
    simData <- simOneTheta(n = 20000,  # 1000000 in original script
                           tree = tree, colSpecies = colSpecies,
                           theta = params["theta"], sigma2 = params["sigma2"],
                           alpha = params["alpha"], beta = params["beta"])

    # run the evemodel to test for regulatory activity shift on the selected branch(es)
    fitEVEsim <- twoThetaTest(tree = tree, gene.data = simData,
                              isTheta2edge = isTheta2edge, colSpecies = colSpecies,
                              cores = 6)

    # Compute empirical p-values from simulation and write the TwoThetaTest results to file
    shift <- fitEVE$twoThetaRes$par[,"theta2"] - fitEVE$twoThetaRes$par[,"theta1"]
    evePvalues <- qvalue::empPvals(fitEVE$LRT, fitEVEsim$LRT)
    resulteve <- data.frame(gene.id = rownames(data_norm),
                            tissue = tissue,
                            branch = branch,
                            pvalue = evePvalues,
                            adjpvalue = p.adjust(evePvalues, method = 'BH'),
                            logFC = shift,
                            baseMean = fitEVE$twoThetaRes$par[,"theta1"],
                            method = "EVEemp")

    all_results <- rbind(all_results, resulteve)
    write.table(all_results,
                here(results_directory, paste0(data_name, ".tsv")),
                sep='\t', quote = FALSE)
    
    if (tissue == "Liver") {
      ################################################################################
      ## Design - batch effect
      ################################################################################
      berthelot <- grepl("berthelot", colnames(data_norm)) + 0
      design_formula_b <- as.formula(~ berthelot + condition)
      design_data_b <- data.frame(condition = condition, berthelot = berthelot)
      rownames(design_data_b) <- colnames(data_norm)
      design_data_b$condition <- factor(design_data_b$condition)
      design_data_b$berthelot <- factor(design_data_b$berthelot)
      design_batch <- model.matrix(design_formula_b, design_data_b)

      ################################################################################
      ## phyloDE - BM - batch effect
      ################################################################################
      fitlimmaBM <- phyloDE::phylolmFit(data_norm,
                                           design = design_batch,
                                           phy = tree,
                                           model = 'BM',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE,
                                           ddf_method = "Samples",
                                           medianOU = FALSE,
                                           ncores = 6)

      fitlimmaBM <- limma::eBayes(fitlimmaBM, trend = TRUE)
      limmaBMPvalues <- fitlimmaBM$p.value[, ncol(fitlimmaBM$p.value)]

      resultlimmaBM <- data.frame(gene.id = rownames(data_norm),
                                  tissue = tissue,
                                  branch = branch,
                                  pvalue = limmaBMPvalues,
                                  adjpvalue = p.adjust(limmaBMPvalues, method = 'BH'),
                                  logFC = fitlimmaBM$coefficients[, ncol(fitlimmaBM$coefficients)],
                                  baseMean = fitlimmaBM$coefficients[, 1],
                                  method = "phyloDEBMbatch")

      all_results_batch <- rbind(all_results_batch, resultlimmaBM)
      write.table(all_results_batch,
                  here(results_directory, paste0(data_name, "_batch.tsv")),
                  sep='\t', quote = FALSE)

      ################################################################################
      ## phyloDE - OU
      ################################################################################
      fitlimmaOU <- phyloDE::phylolmFit(data_norm,
                                           design = design_batch,
                                           phy = tree,
                                           model = 'OUfixedRoot',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE,
                                           ddf_method = "Samples",
                                           medianOU = FALSE,
                                           ncores = 6)

      fitlimmaOU <- limma::eBayes(fitlimmaOU, trend = TRUE)
      limmaOUPvalues <- fitlimmaOU$p.value[, ncol(fitlimmaOU$p.value)]

      resultlimmaOU <- data.frame(gene.id = rownames(data_norm),
                                  tissue = tissue,
                                  branch = branch,
                                  pvalue = limmaOUPvalues,
                                  adjpvalue = p.adjust(limmaOUPvalues, method = 'BH'),
                                  logFC = fitlimmaOU$coefficients[, ncol(fitlimmaOU$coefficients)],
                                  baseMean = fitlimmaOU$coefficients[, 1],
                                  method = "phyloDEOUbatch")

      all_results_batch <- rbind(all_results_batch, resultlimmaOU)
      write.table(all_results_batch,
                  here(results_directory, paste0(data_name, "_batch.tsv")),
                  sep='\t', quote = FALSE)
    } else if (tissue == "Heart") {
      ################################################################################
      ## Design - batch effect
      ################################################################################
      faulkes <- grepl("faulkes", colnames(data_norm)) + 0
      bens <- grepl("bens", colnames(data_norm)) + 0
      design_formula_b <- as.formula(~ bens + faulkes + condition)
      design_data_b <- data.frame(condition = condition, faulkes = faulkes, bens = bens)
      rownames(design_data_b) <- colnames(data_norm)
      design_data_b$condition <- factor(design_data_b$condition)
      design_data_b$faulkes <- factor(design_data_b$faulkes)
      design_data_b$bens <- factor(design_data_b$bens)
      design_batch <- model.matrix(design_formula_b, design_data_b)

      ################################################################################
      ## phyloDE - BM - batch effect
      ################################################################################
      fitlimmaBM <- phyloDE::phylolmFit(data_norm,
                                           design = design_batch,
                                           phy = tree,
                                           model = 'BM',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE,
                                           ddf_method = "Samples",
                                           medianOU = FALSE,
                                           ncores = 6)

      fitlimmaBM <- limma::eBayes(fitlimmaBM, trend = TRUE)
      limmaBMPvalues <- fitlimmaBM$p.value[, ncol(fitlimmaBM$p.value)]

      resultlimmaBM <- data.frame(gene.id = rownames(data_norm),
                                  tissue = tissue,
                                  branch = branch,
                                  pvalue = limmaBMPvalues,
                                  adjpvalue = p.adjust(limmaBMPvalues, method = 'BH'),
                                  logFC = fitlimmaBM$coefficients[, ncol(fitlimmaBM$coefficients)],
                                  baseMean = fitlimmaBM$coefficients[, 1],
                                  method = "phyloDEBMbatch")

      all_results_batch <- rbind(all_results_batch, resultlimmaBM)
      write.table(all_results_batch,
                  here(results_directory, paste0(data_name, "_batch.tsv")),
                  sep='\t', quote = FALSE)

      ################################################################################
      ## phyloDE - OU
      ################################################################################
      fitlimmaOU <- phyloDE::phylolmFit(data_norm,
                                           design = design_batch,
                                           phy = tree,
                                           model = 'OUfixedRoot',
                                           measurement_error = TRUE,
                                           use_consensus = TRUE,
                                           ddf_method = "Samples",
                                           medianOU = FALSE,
                                           ncores = 6)

      fitlimmaOU <- limma::eBayes(fitlimmaOU, trend = TRUE)
      limmaOUPvalues <- fitlimmaOU$p.value[, ncol(fitlimmaOU$p.value)]

      # go <- limma::goana(fitlimmaOU, coef=2, species = "Mm",
      #                    geneid = sub("ENSMUSG[0]*", "", rownames(fitlimmaOU)))
      # topGO(go, p.value = 0.1)
      # 
      # kk <- limma::kegga(fitlimmaOU, coef=2, species = "Mm",
      #                    geneid = sub("ENSMUSG[0]*", "", rownames(fitlimmaOU)))
      # topKEGG(kk, p.value = 0.1)
      
      aa <- sapply(Lkeys(org.Mm.eg.db::org.Mm.egGO2ALLEGS)[1:100], function(x) grep(x, rownames(fitlimmaOU)))
      aa <- do.call(c, aa)
      resultlimmaOU <- data.frame(gene.id = rownames(data_norm),
                                  tissue = tissue,
                                  branch = branch,
                                  pvalue = limmaOUPvalues,
                                  adjpvalue = p.adjust(limmaOUPvalues, method = 'BH'),
                                  logFC = fitlimmaOU$coefficients[, ncol(fitlimmaOU$coefficients)],
                                  baseMean = fitlimmaOU$coefficients[, 1],
                                  method = "phyloDEOUbatch")

      all_results_batch <- rbind(all_results_batch, resultlimmaOU)
      write.table(all_results_batch,
                  here(results_directory, paste0(data_name, "_batch.tsv")),
                  sep='\t', quote = FALSE)
    }
  }
}

################################################################################
## COMPARISON : EVE tree myr vs EVE tree subs
################################################################################
tissue <- "Liver" # one of "Liver" or "Heart"
branch <- "hgla"  # one of "MR", "hgla", "fdam" or "cpor"
Branch <- "Hgla"  # one of "MR", "hgla", "fdam" or "cpor"

## EVE with tree Hassler (this script)
all_results <- read.table(here(results_directory, paste0(data_name, ".tsv")))
results_eve <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "EVEemp",]
results_eve <- results_eve[complete.cases(results_eve), ]
dim(results_eve)
## EVE with tree Parey (reproduce script)
results_directory_reproduce <- here(paste0(data_file, "_analysis_reproduce_", "2025-10-29"))
resulteve <- read.table(here(results_directory_reproduce, paste0(data_name, "_", tissue, "_", branch, "_EVE.tsv")))
resulteve <- resulteve[complete.cases(resulteve), ]
dim(resulteve)

res_both <- merge(results_eve, resulteve, by = "gene.id")
hist(res_both$empirical.pval - res_both$pvalue)

sel_paper <- resulteve[resulteve$empirical.pval <= 0.05, "gene.id"]
sel_eve <- results_eve[results_eve$pvalue <= 0.05, "gene.id"]
length(sel_paper)
length(sel_eve)
# all genes in paper are recovered in eve analysis
anyNA(match(sel_paper, sel_eve))
# genes in new eve analysis but not in paper
results_eve[match(sel_eve[which(is.na(match(sel_eve, sel_paper)))], results_eve$gene.id), ]
# same genes are selected


################################################################################
## COMPARISON : EVE vs phyloDE
################################################################################
tissue <- "Heart" # one of "Liver" or "Heart"
branch <- "MR"  # one of "MR", "hgla", "fdam" or "cpor"

all_results <- read.table(here(results_directory, paste0(data_name, ".tsv")))
all_results_batch <- read.table(here(results_directory, paste0(data_name, "_batch.tsv")))
all_results <- rbind(all_results, all_results_batch)
## EVE
results_eve <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "EVEemp",]
results_eve <- results_eve[complete.cases(results_eve), ]
dim(results_eve)
hist(results_eve$pvalue)
## BM
results_bm <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "phyloDEBM",]
results_bm <- results_bm[complete.cases(results_bm), ]
dim(results_bm)
hist(results_bm$pvalue)
## OU
results_ou <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "phyloDEOU",]
results_ou <- results_ou[complete.cases(results_ou), ]
dim(results_ou)
hist(results_ou$pvalue)
## BM batch
results_bm_batch <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "phyloDEBMbatch",]
results_bm_batch <- results_bm_batch[complete.cases(results_bm_batch), ]
dim(results_bm_batch)
hist(results_bm_batch$pvalue)
## OU
results_ou_batch <- all_results[all_results$tissue == tissue & all_results$branch == branch & all_results$method == "phyloDEOUbatch",]
results_ou_batch <- results_ou_batch[complete.cases(results_ou_batch), ]
dim(results_ou_batch)
hist(results_ou_batch$pvalue)
sum(results_ou_batch$adjpvalue <= 0.05)

## upsetr
library(UpSetR)
x = list (
  eve = results_eve[results_eve$pvalue <= 0.05, "gene.id"],
  eveadj = results_eve[results_eve$adjpvalue <= 0.05, "gene.id"],
  # bm = results_bm[results_bm$adjpvalue <= 0.05, "gene.id"],
  # ou = results_ou[results_ou$adjpvalue <= 0.05, "gene.id"],
  # bmbatch = results_bm_batch[results_bm_batch$adjpvalue <= 0.05, "gene.id"],
  oubatch = results_ou_batch[results_ou_batch$adjpvalue <= 0.05, "gene.id"]
)
upset(fromList(x), nsets = 4, order.by = "freq")
ggsave(filename = here(results_directory, paste0(data_name, "_upset.pdf")))

results_ou_batch[results_ou_batch$adjpvalue <= 0.05, ]
