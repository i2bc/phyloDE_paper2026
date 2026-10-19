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
# using the mammals dataset
################################################################################

## INPUT #######################################################################
# Formatted data
# Saved in `mammals_RNAseq/results`
################################################################################

## OUTPUT ######################################################################
# Simulated data
# Saved in `mammals_RNAseq/results/YYY-MM-DD_simulations_mammals_RNAseq`
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
data_chen   <- readRDS(here("mammals_RNAseq", "results", "chen2019_primates_cpd.rds"))
counts_noNA <- data_chen@count.matrix
length_noNA <- data_chen@length.matrix
colData     <- data_chen@sample.annotations
colnames(colData)[colnames(colData) == "condition"] <- "primates"
colData$primates[colData$species == "marmoset"] <- 2 # marmoset miss annotated ?
colData$rodents <- rep(1, nrow(colData))
colData$rodents[grep("mus|rat", colData$species)] <- 2
tree        <- data_chen@tree

################################################################################
## Files management
################################################################################

# here_dir <- sub("save", "work", here())
here_dir        <- here("mammals_RNAseq")
datestamp_day   <- format(Sys.time(), "%Y-%m-%d")
simus_directory <- paste0(datestamp_day, "_simulations_mammals")
simus_directory <- file.path(here_dir, "results", simus_directory)
dir.create(simus_directory, showWarnings = FALSE)

################################################################################
## run DESeq2 to get parameters estimation
################################################################################

# Correct by families
dds <- DESeqDataSetFromMatrix(counts_noNA, colData, ~ primates)

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

dispersions_chen <- DESeq2::dispersions(dds)
seqdepth_chen    <- mean(colSums(counts_noNA))
minfact_chen     <- min(colSums(counts_noNA)/median(colSums(counts_noNA)))
maxfact_chen     <- max(colSums(counts_noNA)/median(colSums(counts_noNA)))
relmeans_chen_p  <- rowMeans(counts_noNA[,colData$primates == 2])
relmeans_chen_r  <- rowMeans(counts_noNA[,colData$rodents  == 2])
n.vars           <- dim(counts_noNA)[1]

to_remove        <- is.na(dispersions_chen)
n.vars           <- n.vars - sum(to_remove) 
dispersions_chen <- dispersions_chen[!to_remove]
relmeans_chen     <- list(primates = relmeans_chen_p[!to_remove],
                          rodents  = relmeans_chen_r[!to_remove])
plot(log(relmeans_chen$rodents), log(dispersions_chen))


################################################################################
## Format Tree
################################################################################

ape::is.ultrametric(tree) # check ultrametric

N <- length(tree$tip.label) # total number of observations (with replicates)
tree$edge.length <- tree$edge.length / max(diag(vcv(tree))[1:N]) # normalize tree to unit height

id_species <- sub("\\..*", "", tree$tip.label)
id_species <- factor(id_species)
names(id_species) <- tree$tip.label

################################################################################
## Conditions design
################################################################################

all_cond_types <- c("primates", "rodents")
all_conds <- sapply(all_cond_types, function(x) {
  setNames(colData[[x]], rownames(colData))
},USE.NAMES=TRUE, simplify=FALSE)

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

params$dataset <- paste("mammals", params$condType, params$modelProcess, 
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
    seqdepth                   = seqdepth_chen,
    minfact                    = minfact_chen,
    maxfact                    = maxfact_chen,
    relmeans                   = relmeans_chen[[cond_type]],
    dispersions                = dispersions_chen * fact_disp,
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
## Save transformed tree
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

tree_tr <- transform_tree_model_OUfixedRoot(tree, 0.7, log(2)/0.5)
nEff    <- list(primates = nEffNaive(tree_tr, all_conds[["primates"]]), 
                rodents  = nEffNaive(tree_tr, all_conds[["rodents"]]))
saveRDS(list(treeTr = tree_tr, nEff = nEff, allConds = all_conds), file.path(simus_directory, "designs.rds"))