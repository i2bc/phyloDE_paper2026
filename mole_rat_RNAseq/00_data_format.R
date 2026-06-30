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
# Program to format the data from:
# Daunesse, M., Liorzou, E., Parey, E., Villar, D. and Berthelot, C. (2025), 
# 'Coordinated shifts in gene expression and regulation during mole-rat evolution', 
# bioRxiv, 2025.11.24.690111.
################################################################################

## INPUT #######################################################################
# Data from https://gitlab.com/maelle.daunesse/daunesse_et_al_2025
# Saved in `mole_rat_RNAseq/data/Daunesse2025`
################################################################################

## OUTPUT ######################################################################
# Formatted data as a compcodeR object
# Saved in `mole_rat_RNAseq/data/data_formatted`
################################################################################


################################################################################
################################################################################

library(here)
library(compcodeR)
library(phyloDE)

################################################################################
## Parameters
################################################################################
data_file <- "mole_rat_RNAseq/data"
tissue <- "Liver" # one of "Liver" or "Heart"
data_name <- paste0("mole_rat_RNAseq_data_", tissue)
dir.create(here(data_file, "data_formatted"), showWarnings = FALSE)

################################################################################
## Read and format count and length data
################################################################################
# Read the raw counts
raw_counts <- read.table(here(data_file, "Daunesse2025", paste0(tissue, "_salmon_matrix_raw_counts.tab")))
raw_counts <- round(raw_counts) # expected RSEM values -> counts

# Read length information
leng <- read.table(here(data_file, "Daunesse2025", "lengthMat.txt"))

# remove unknown lengths
length_noNA <- leng[complete.cases(leng),]
raw_counts_noNA <- raw_counts[match(rownames(length_noNA), rownames(raw_counts)), ]
table(rownames(raw_counts_noNA) == rownames(length_noNA))

# find species for each sample
species_vector <- sapply(strsplit(colnames(raw_counts_noNA), "_"), function(x) x[1])

# replicates columns in lengths
length_noNA <- length_noNA[, species_vector]
colnames(length_noNA) <- colnames(raw_counts_noNA)

################################################################################
## Read and format condition data
################################################################################

# Design phyloDE
conditionMR <- grepl("Hgla|Fdam", colnames(raw_counts_noNA)) + 1
conditionHgla <- grepl("Hgla", colnames(raw_counts_noNA)) + 1
conditionFdam <- grepl("Fdam", colnames(raw_counts_noNA)) + 1
conditionCpor <- grepl("Cpor", colnames(raw_counts_noNA)) + 1

colData <- data.frame(condition = as.factor(conditionMR),
                      conditionMR = as.factor(conditionMR),
                      conditionHgla = as.factor(conditionHgla),
                      conditionFdam = as.factor(conditionFdam),
                      conditionCpor = as.factor(conditionCpor),
                      species = species_vector,
                      id = colnames(raw_counts_noNA))
rownames(colData) <- colData$id

################################################################################
## Read and format Tree
################################################################################
library(ape)
# get ultrametric mammal tree Hassler et al. 2020 (https://doi.org/10.1080/01621459.2020.1799812).
mammals_tree_ultra <- read.tree(file = url("https://raw.githubusercontent.com/suchard-group/incomplete_measurements/master/data/mammals_newick.txt"))
# Trim to mole rats
species_parey <- c("Mmus", "Mus_musculus",
                   "Hgla", "Heterocephalus_glaber",
                   "Fdam", "Cryptomys_damarensis",
                   "Cpor", "Cavia_porcellus")
species_parey <- matrix(species_parey, ncol = 2, byrow = TRUE)
colnames(species_parey) <- c("common", "scientific")
mammals_tree_ultra_trimmed <- keep.tip(mammals_tree_ultra, species_parey[, "scientific"])
mammals_tree_ultra_trimmed$tip.label <- species_parey[match(mammals_tree_ultra_trimmed$tip.label, species_parey[, "scientific"]), "common"]
# mammals_tree_ultra_trimmed$edge.length <- mammals_tree_ultra_trimmed$edge.length/mean(diag(vcv(mammals_tree_ultra_trimmed)))
mammals_tree_ultra_trimmed <- ladderize(mammals_tree_ultra_trimmed, right = TRUE)
write.tree(mammals_tree_ultra_trimmed, file = here("mole_rat_RNAseq/data/data_formatted/species_tree_mya.nwk"))
plot(mammals_tree_ultra_trimmed)
edgelabels(round(mammals_tree_ultra_trimmed$edge.length, 2))

## Compare with Parey tree
treeParey <- read.tree(file = here("mole_rat_RNAseq/data/Daunesse2025/species_tree.nwk"))
# treeParey$edge.length <- treeParey$edge.length/mean(diag(vcv(treeParey)))
plot(treeParey)
edgelabels(round(treeParey$edge.length, 2))

## choose tree
tree <- mammals_tree_ultra_trimmed

## Species names
# Match
tree_data_cor <- match(tree$tip.label, colData$species)
data_tree_cor <- match(colData$species, tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
colData$species[is.na(data_tree_cor)]

## Format Tree
# Get rid of species not in data
tree <- drop.tip(tree, tip = tree$tip.label[is.na(tree_data_cor)])
plot(tree)
# add replicates
tree_rep <- phyloDE::addReplicatesOnTree(tree, colData)
plot(tree_rep)

## Reorder data
corr <- match(tree_rep$tip.label, colData$id)
colData <- colData[corr, ]
raw_counts_noNA <- raw_counts_noNA[, corr]
length_noNA <- length_noNA[, corr]

################################################################################
# Format for compcodeR
################################################################################
## id species 
colData$id.species <- colData$species
## create object
cpd <- phyloCompData(tree = tree_rep,
                     count.matrix = raw_counts_noNA, 
                     sample.annotations = colData, 
                     info.parameters = list(dataset = data_name, uID = "1"),
                     length.matrix = as.matrix(length_noNA))
## check obj conformity
check_phyloCompData(cpd)

## save data to rds file
saveRDS(cpd, file = here(data_file, "data_formatted", data_name))
