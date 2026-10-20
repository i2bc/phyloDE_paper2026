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
################################################################################

## OUTPUT ######################################################################
# Raw count and length data 
# Saved in `mammals_RNAseq/data/`
################################################################################

library(here)

## read file of existing TPM exp matrix for liver 
tab <- read.csv2("https://raw.githubusercontent.com/evee-tools/EVEE-Tools/master/data/liver.exp_matrix.txt",sep="\t")
head(tab)
str(tab)
rownames(tab) <- tab$X.OG
tabRaw <- lengthRaw <- effLengthRaw <- TPMRaw <- FPKMRaw <- tab
tabRaw[,] <- lengthRaw[,] <- effLengthRaw[,] <- TPMRaw[,] <- FPKMRaw[,]  <- NA

## read file associating species to file directory containing exp data
rsem.list <- read.table("https://raw.githubusercontent.com/evee-tools/EVEE-Tools/master/data/liver.rsem.list",sep="\t")
head(rsem.list)
colnames(rsem.list) <- c("species","species.id","path")
## add unique id for sample
rownames(rsem.list) <- colnames(tab)[-c(1,2)]
head(rsem.list)

## read orthogroups 
orthogroups <- read.csv2("https://raw.githubusercontent.com/evee-tools/EVEE-Tools/master/data/mammals.ortholog_groups.txt",
                         sep = "\t", header = FALSE, row.names = 1)

## Process data
for(samid in rownames(rsem.list)[1:2]){
  print(samid)
  ## read expression data for current sample
  curpath <- rsem.list[samid, "path"]
  curspeciesid <- rsem.list[samid, "species.id"]
  rsem.file <- read.table(paste0("https://raw.githubusercontent.com/evee-tools/EVEE-Tools/master/", curpath), header = TRUE)
  i <- 0
  ## get expression for each OG
  for(idOG in rownames(tab)){
    i <- i+1
    if(i%%1000==0){print(i)}
    ## find OG in species orthogroup
    idg <- grep(curspeciesid, orthogroups[idOG, ])
    if(length(idg) > 0){ 
      GN <- strsplit(orthogroups[idOG, idg], "[|]")[[1]][2]
      selin <- rsem.file$gene_id %in% GN
      if(length(table(selin)) > 1){
        TPMRaw[idOG, samid] <- rsem.file[selin,"TPM"] 
        FPKMRaw[idOG, samid] <- rsem.file[selin,"FPKM"] 
        tabRaw[idOG, samid] <- rsem.file[selin,"expected_count"]
        lengthRaw[idOG, samid] <- rsem.file[selin,"length"]
        effLengthRaw[idOG, samid] <- rsem.file[selin,"effective_length"]
      }
    }
  }
}

## update OG name and gene name
tabRaw$X.OG <- lengthRaw$X.OG <- effLengthRaw$X.OG <- rownames(tabRaw)
TPMRaw$X.OG <- FPKMRaw$X.OG <- rownames(tabRaw)
tabRaw$GENE_NAME <- lengthRaw$GENE_NAME <- effLengthRaw$GENE_NAME <- tab$GENE_NAME
TPMRaw$GENE_NAME <- FPKMRaw$GENE_NAME <- tab$GENE_NAME

## write file
write.table(tabRaw, file = here("mammals_RNAseq/data/rawCounts.txt"))
write.table(lengthRaw, file = here("mammals_RNAseq/data/rawLengths.txt"))
write.table(effLengthRaw, file = here("mammals_RNAseq/data/rawEffLengths.txt"))
write.table(TPMRaw, file = here("mammals_RNAseq/data/rawTPM.txt"))
write.table(FPKMRaw, file = here("mammals_RNAseq/data/rawFPKM.txt"))


