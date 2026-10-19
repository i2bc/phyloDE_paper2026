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
# Program to simulate data from empirical parameters
# using the crayfish dataset
################################################################################

## INPUT #######################################################################
# Formatted data
# Saved in `crayfish_RNAseq/results`
################################################################################

## OUTPUT ######################################################################
# Simulated data
# Saved in `crayfish_RNAseq/results/YYY-MM-DD_simulations_crayfish_RNAseq`
################################################################################

library(here)
library(compcodeR)
library(ape)
library(DESeq2)
library(parallel)

NB_THREADS <- 10

################################################################################
## Load Data
################################################################################

## read data
data_stern  <- readRDS(here("crayfish_RNAseq", "results", "stern2018_cpd.rds"))
counts_noNA <- data_stern@count.matrix
length_noNA <- data_stern@length.matrix
colData     <- data_stern@sample.annotations
colnames(colData)[colnames(colData) == "condition"] <- "sights"
tree        <- data_stern@tree

################################################################################
## Files management
################################################################################

# here_dir <- sub("save", "work", here())
here_dir        <- here("crayfish_RNAseq")
datestamp_day   <- format(Sys.time(), "%Y-%m-%d")
simus_directory <- paste0(datestamp_day, "_simulations_crayfish")
simus_directory <- file.path(here_dir, "results", simus_directory)
dir.create(simus_directory, showWarnings = FALSE)

################################################################################
## run DESeq2 to get parameters estimation
################################################################################

# Correct by families
dds <- DESeqDataSetFromMatrix(counts_noNA, colData, ~ sights)

# normalisation factors for library size
dds      <- estimateSizeFactors(dds)
size_fac <- sizeFactors(dds)

# compute normalization with length
mat_size_fac <- matrix(size_fac, ncol=length(size_fac), nrow=length(counts_noNA[,1]), byrow=T)
normFactors  <- (mat_size_fac*length_noNA) / exp(rowMeans(log(mat_size_fac*length_noNA)))
normalizationFactors(dds) <- as.matrix(normFactors)

# analysis
dds <- DESeq2::DESeq(dds)
dds <- DESeq2::replaceOutliersWithTrimmedMean(dds)
dds <- DESeq2::DESeq(dds)

################################################################################
## Get empirical moments for simulation
################################################################################

dispersions_stern <- DESeq2::dispersions(dds)
seqdepth_stern    <- mean(colSums(counts_noNA))
minfact_stern     <- min(colSums(counts_noNA)/median(colSums(counts_noNA)))
maxfact_stern     <- max(colSums(counts_noNA)/median(colSums(counts_noNA)))
relmeans_stern    <- rowMeans(counts_noNA[,colData$sights == 2])
n.vars            <- dim(counts_noNA)[1]

anyNA(dispersions_stern) # no NA -> continue
plot(log(relmeans_stern), log(dispersions_stern))

################################################################################
## Format Tree
################################################################################

ape::is.ultrametric(tree) # check ultrametric

n <- length(unique(sub("_.*", "", tree$tip.label))) # number of species
N <- length(tree$tip.label) # total number of observations (with replicates)

tree$edge.length <- tree$edge.length / max(diag(vcv(tree))[1:N]) # normalize tree to unit height

id_species <- sub("_.*", "", tree$tip.label)
id_species <- factor(id_species)
names(id_species) <- tree$tip.label

species_names <- sapply(1:length(colnames(length_noNA)), function(i) strsplit(colnames(length_noNA),"_")[[i]][1])
species_names <- sapply(1:length(species_names), function(i) strsplit(species_names,"[[:upper:]]$")[[i]][1])

################################################################################
## Conditions design
################################################################################

all_cond_species <- NULL
all_cond_types   <- NULL
all_conds        <- NULL

## Blocs
species_names         <- unique(sub("_.*", "", tree$tip.label))
cond_species          <- rep(c(1, 2), each = length(species_names) / 2)
names(cond_species)   <- species_names
cond_species["CHAMU"] <- 1
id_cond               <- id_species
id_cond               <- cond_species[as.vector(id_cond)]
id_cond               <- as.factor(id_cond)
names(id_cond)        <- tree$tip.label

plot(tree)
tiplabels(pch = 21, col = id_cond, bg = id_cond)

all_cond_types              <- c(all_cond_types, "bloc")
all_conds[["bloc"]]         <- id_cond
all_cond_species[["bloc"]] <- cond_species

################################################################################
## Parameters for the simulation
################################################################################

params <- expand.grid(
  rep               = seq_len(20),
  effectSize        = 3,
  selectionStrength = log(2) / 0.5,
  propVarTree       = 0.7,
  factDisp          = 1,
  modelProcess      = c("OU"),
  condType          = all_cond_types,
  stringsAsFactors  = FALSE
)

params$dataset <- paste("crayfish", params$condType, params$modelProcess, 
                        params$effectSize, params$selectionStrength, 
                        params$propVarTree, params$factDisp, sep = "_") 

saveRDS(params, file.path(simus_directory, "params.rds"))

################################################################################
## Simulations
################################################################################

invisible(mclapply(1:nrow(params), function(idx) {
  
  rep                <- params$rep[idx]
  effect_size        <- params$effectSize[idx]
  selection_strength <- params$selectionStrength[idx]
  prop_var_tree      <- params$propVarTree[idx]
  fact_disp          <- params$factDisp[idx]
  model_process      <- params$modelProcess[idx]
  cond_type          <- params$condType[idx]
  dataset            <- params$dataset[idx]
  
  set.seed(18570823 + idx)
  
  simus_obj <- generateSyntheticData(
    dataset                    = dataset,
    n.vars                     = n.vars,
    samples.per.cond           = N / 2,
    n.diffexp                  = 250,
    repl.id                    = rep,
    seqdepth                   = seqdepth_stern,
    minfact                    = minfact_stern,
    maxfact                    = maxfact_stern,
    relmeans                   = relmeans_stern,
    dispersions                = dispersions_stern * fact_disp,
    fraction.upregulated       = 0.5,
    between.group.diffdisp     = FALSE,
    filter.threshold.total     = 1,
    filter.threshold.mediancpm = 0,
    fraction.non.overdispersed = 0,
    random.outlier.high.prob   = 0,
    random.outlier.low.prob    = 0,
    single.outlier.high.prob   = 0,
    single.outlier.low.prob    = 0,
    effect.size                = effect_size,
    output.file                = file.path(simus_directory, paste0(dataset, "_", rep, ".rds")),
    tree                       = tree,
    prop.var.tree              = prop_var_tree,
    id.condition               = all_conds[[cond_type]],
    id.species                 = id_species,
    model.process              = model_process,
    selection.strength         = ifelse(model_process == "OU", selection_strength, 0)
  )
  
}, mc.cores = NB_THREADS))

################################################################################
## Select 2*k species with maximum nEff
################################################################################

# Rescale tree to unit height
rescale_tree <- function(tree) {
  tree$edge.length <- tree$edge.length / tree_height(tree)
  return(tree)
}

# Get maximum tree height
tree_height <- function(tree) {
  return(max(ape::node.depth.edgelength(tree)))
}

# Transform tree under Brownian Motion model with errors
transform_tree_model_BM <- function(tree, lambda) {
  return(rescale_tree(phylolm::transf.branch.lengths(tree, "lambda", parameters = list(lambda = lambda))$tree))
}

# Transform tree under Ornstein-Uhlenbeck model with errors
transform_tree_model_OUfixedRoot <- function(tree, lambda, alpha) {
  tree <- phylolm::transf.branch.lengths(tree, "OUfixedRoot", parameters = list(alpha = alpha))$tree
  return(rescale_tree(phylolm::transf.branch.lengths(tree, "lambda", parameters = list(lambda = lambda))$tree))
}

# Calculate effective sample size under phylogenetic model
# Adapted from compcodeR package
nEffNaive <- function(tree, id.condition) {
  if (length(unique(id.condition)) != 2) stop("There should be only two conditions")
  if (is.null(tree)) {
    n           <- length(id.condition)
    id_cond_uni <- id.condition
    Vinv        <- diag(1, n)
  } else {
    id.condition <- checkParamVector(id.condition, "id.condition", tree)
    tree_uni     <- tree
    id_cond_uni  <- id.condition[tree_uni$tip.label]
    n            <- length(tree_uni$tip.label)
    V            <- ape::vcv(tree_uni)
    Vinv         <- try(solve(V))
    if (inherits(Vinv, 'try-error')) Vinv <- ginv(V)
  }
  # Intercept
  X <- rep(1, n)
  # Model matrix
  Z <- stats::model.matrix(~factor(id_cond_uni))[, -1, drop = FALSE]
  # Variance of the condition
  X <- cbind(X, Z)
  nEff <- 1 / solve(t(X) %*% Vinv %*% X)[2, 2]
  return(nEff)
}

# Check and reorder condition vector to match tree tip labels
# Adapted from compcodeR package
checkParamVector <- function(x, name, tree) {
  N <- length(tree$tip.label)
  
  if (length(x) != N) {
    stop(paste0("`", name, "` should have the same length as the number of taxa in the tree."))
  }
  if (is.null(tree$tip.label)) stop("The tips of the phylogeny are not named.")
  if (is.null(names(x))){
    warning(paste0("`", name, "` is not named. I'm naming them, assuming they are in the same order as the tree."))
    names(x) <- tree$tip.label
  } else {
    if (!all(tree$tip.label == names(x))){
      correspondances <- match(tree$tip.label, names(x))
      if (length(na.omit(unique(correspondances))) != length(tree$tip.label)){
        stop(paste0("`", name, "` names do not match the tip labels."))
      }
      warning(paste0("`", name, "` was not sorted in the correct order, when compared with the tips label. I am re-ordering it."))
      x <- x[correspondances]
    }
  }
  return(x)
}

#Extract subtree containing only specified species (with all their replicates) 
subset_tree_by_species <- function(tree, species_names) {
  
  tip_species  <- sub("_.*", "", tree$tip.label)
  tips_to_keep <- tree$tip.label[tip_species %in% species_names]
  
  if (length(tips_to_keep) == 0) {
    stop("No tips found matching the provided species names")
  }
  subtree <- ape::keep.tip(tree, tips_to_keep)
  
  return(rescale_tree(subtree))
}

# This function finds the optimal subset of k species per condition that
# maximizes the effective sample size (nEff) under a phylogenetic model.
# The function performs an exhaustive search over all possible combinations
# of k species in condition 1 and k species in condition 2.
select_optimal_species <- function(tree, cond_species, k, model = "BM",
                                   lambda = 0.6, alpha = log(2) / 0.5) {
  
  species_cond1 <- names(cond_species)[cond_species == 1]
  species_cond2 <- names(cond_species)[cond_species == 2]
  
  if (k > length(species_cond1)) {
    stop(paste0("k (", k, ") is larger than the number of species in condition 1 (",
                length(species_cond1), ")"))
  }
  if (k > length(species_cond2)) {
    stop(paste0("k (", k, ") is larger than the number of species in condition 2 (",
                length(species_cond2), ")"))
  }
  
  comb_cond1 <- combn(species_cond1, k, simplify = FALSE)
  comb_cond2 <- combn(species_cond2, k, simplify = FALSE)
  
  best_neff      <- -Inf
  best_species_1 <- NULL
  best_species_2 <- NULL
  
  # Test all combinations
  for (sp1 in comb_cond1) {
    for (sp2 in comb_cond2) {
      
      selected_species   <- c(sp1, sp2)
      subtree            <- subset_tree_by_species(tree, selected_species)
      subtree_cond       <- cond_species[selected_species]
      tip_species        <- sub("_.*", "", subtree$tip.label)
      id_cond_sub        <- subtree_cond[tip_species]
      names(id_cond_sub) <- subtree$tip.label
      
      if (model == "BM") {
        tree_transformed <- transform_tree_model_BM(subtree, lambda)
      } else {
        tree_transformed <- transform_tree_model_OUfixedRoot(subtree, lambda, alpha)
      }
      
      neff   <- nEffNaive(tree_transformed, id_cond_sub)
      
      if (neff > best_neff) {
        best_neff      <- neff
        best_species_1 <- sp1
        best_species_2 <- sp2
        best_tree      <- subtree
        best_tree_tr   <- tree_transformed
      }
    }
  }
  
  return(list(
    species_cond1 = best_species_1,
    species_cond2 = best_species_2,
    neff          = best_neff,
    tree          = best_tree,
    treeTr        = best_tree_tr
  ))
} 

# For each combination of parameters (selection strength, phylogenetic signal,
# evolutionary model, condition type) and each value of k from 1 to k_max,
# find the species subset (k per condition) that maximizes nEff
design_params <- unique(params[,colnames(params) %in% c("selectionStrength", "propVarTree", "modelProcess", "condType")])
kmax          <- 6

designs <- mclapply(1:nrow(design_params), function(idx) {
  lapply(1:kmax, function(k) {
    select_optimal_species(
      tree,
      all_cond_species[[design_params$condType[idx]]],
      k,
      design_params$modelProcess[idx],
      design_params$propVarTree[idx],
      design_params$selectionStrength[idx]
    )
  })
}, mc.cores = NB_THREADS)

saveRDS(list(
  designs = designs,
  params  = design_params,
  kmax    = kmax
), file.path(simus_directory, "designs.rds"))