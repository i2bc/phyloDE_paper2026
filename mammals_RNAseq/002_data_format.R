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
# Chen, J., Swofford, R., Johnson, J., Cummings, B. B., Rogel, N., Lindblad-Toh, 
# K., Haerty, W., Palma, F. di, & Regev, A. (2019). A quantitative framework for 
# characterizing the evolutionary history of mammalian gene expression. 
# Genome Research, 29(1), 53–63. https://doi.org/10.1101/gr.237636.118
################################################################################

## INPUT #######################################################################
# Data from https://github.com/evee-tools/EVEE-Tools
# Formatted by script `mammals_RNAseq/001_data_format.R`
# Saved in `mole_rat_RNAseq/data`
################################################################################

## OUTPUT ######################################################################
# Formatted data as a compcodeR object
# Saved in `mole_rat_RNAseq/results`
################################################################################


################################################################################
################################################################################

library(here)
library(compcodeR)
library(phyloDE)

################################################################################
## Parameters
################################################################################
data_file <- here("mammals_RNAseq", "data")
result_file <- here("mammals_RNAseq", "results")

################################################################################
## Read and format count and length data
################################################################################
## Read the raw counts
raw_counts_full <- read.table(file.path(data_file, "rawCounts.txt"))
raw_counts <- raw_counts_full[,-c(1,2)]
raw_counts <- round(raw_counts) # expected RSEM values -> counts
## format sample names
vectrep <- paste0(".",sapply(1:dim(raw_counts)[2],function(i) length(strsplit(colnames(raw_counts),"[.]")[[i]]))-1)
colnames(raw_counts) <- paste0(colnames(raw_counts),replace(vectrep,vectrep==".1",""))

## Read length information
leng_full <- read.table(file.path(data_file, "rawLengths.txt"))
leng <- leng_full[,-c(1,2)]
## format sample names
vectrep <- paste0(".",sapply(1:dim(leng)[2],function(i) length(strsplit(colnames(leng),"[.]")[[i]]))-1)
colnames(leng) <- paste0(colnames(leng),replace(vectrep,vectrep==".1",""))

# remove NA
raw_counts_noNA <- raw_counts[complete.cases(raw_counts),]
length_noNA <- leng[complete.cases(raw_counts),colnames(raw_counts_noNA)]

################################################################################
## Read and format condition data
################################################################################

mammalsconds <- read.csv2(file.path(data_file, "speciesgroups.txt"),sep="\t")
mammalsconds
## "primates" : 2= primates, 1= non primates
## "rodents" : 2 = rodents , 1 = non rodents

## for primates
condName <- "primates"
# Species condition (1 = non primates, 2 = primates)
condstp <- mammalsconds$primates
names(condstp) <- mammalsconds$species
# Sample Condition
idSpe <- sapply(1:dim(raw_counts)[2],function(i) strsplit(colnames(raw_counts)[i],"[.]")[[1]][1])
primates <- condstp[idSpe]  
names(primates) <- colnames(raw_counts)
# Format
primates <- factor(primates)
colData <- data.frame(primates)
colData$primates <- factor(colData$primates)
colData$species <- sapply(1:length(rownames(colData)),function(i) strsplit(rownames(colData),"[.]")[[i]][1])
colData$id <- rownames(colData)

################################################################################
## Read and format Tree
################################################################################
library(ape)
# get ultrametric mammal tree Hassler et al. 2020 (https://doi.org/10.1080/01621459.2020.1799812).
mammals_tree_ultra <- read.tree(file = url("https://raw.githubusercontent.com/suchard-group/incomplete_measurements/master/data/mammals_newick.txt"))
species_chen <- c("human", "Homo_sapiens",
                  "chimp", "Pan_troglodytes",
                  "bonobo", "Pan_paniscus",
                  "gorilla", "Gorilla_gorilla",
                  "orangutan", "Pongo_abelii",
                  "rhesus", "Macaca_mulatta",
                  "marmoset", "Callithrix_jacchus",
                  "ferret", "Mustela_putorius",
                  "cow", "Bos_taurus",
                  "musMusculus", "Mus_musculus",
                  "musSpretus", "Mus_spretus",
                  "musCaroli", "Mus_caroli",
                  "rat", "Rattus_norvegicus",
                  "rabbit", "Oryctolagus_cuniculus",
                  "dog", "Canis_lupus",
                  "armadillo", "Dasypus_novemcinctus",
                  "opossum", "Monodelphis_domestica")
species_chen <- matrix(species_chen, ncol = 2, byrow = TRUE)
colnames(species_chen) <- c("common", "scientific")

## Trim
mammals_tree_ultra_trimmed <- keep.tip(mammals_tree_ultra, species_chen[, "scientific"])

## Names
mammals_tree_ultra_trimmed$tip.label <- species_chen[match(mammals_tree_ultra_trimmed$tip.label, species_chen[, "scientific"]), "common"]
mammals_tree_ultra_trimmed <- ladderize(mammals_tree_ultra_trimmed, right = FALSE)

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
## rename key column "condition" in colData
colnames(colData)[grep(condName, colnames(colData))] <- "condition"
## id species 
colData$id.species <- colData$species
## create object
cpd <- phyloCompData(tree = tree_rep,
                     count.matrix = raw_counts_noNA, 
                     sample.annotations = colData, 
                     info.parameters = list(dataset = "chen2019_primates_cpd", uID = "1"),
                     length.matrix = as.matrix(length_noNA))
## check obj conformity
check_phyloCompData(cpd)

## save data to rds file
saveRDS(cpd, file = file.path(result_file, "chen2019_primates_cpd.rds"))
