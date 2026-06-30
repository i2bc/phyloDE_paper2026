library(dplyr)
library(msigdbr)

#Shifts
eve <- read.table("../eve/output/Heart_Hgla_sim1000k/Heart_Hgla_eve_results.csv", header = T, sep = "\t")

# Gene set
gene_sets <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") # category and subcategory are listed in MsigDB

# other possible gene sets
#msigdbr(species = "Mus musculus", category = "C2") # C2
#msigdbr(species = "Mus musculus", category = "H") # Hallmark

# Extracting Up and Down genes
eve_UP <- eve %>%
  filter(empirical.pval < 0.05 & shift > 0)

eve_DO <- eve %>%
  filter(empirical.pval < 0.05 & shift < 0)

# Filtering for genes tested with eve
gene_sets2 <- filter(gene_sets, ensembl_gene %in% eve$gene.id)

# Running hypergeometric test for each gene set
gs_out <- data.frame()
for (GS_name in unique(gene_sets2$gs_name)) {
  gs <- filter(gene_sets2, gs_name == GS_name)
  
  # number of genes tested and number of genes in gene set
  number_genes_tested <- nrow(eve)
  number_gene_in_gs <- nrow(gs)
  
  # Perform hypergeometric test UP
  number_up_genes <- nrow(eve_UP)
  number_up_genes_in_gs <- sum(eve_UP$gene.id %in% gs$ensembl_gene)
  
  Pval <- phyper(q = number_up_genes_in_gs - 1, 
                m = number_up_genes, 
                n = number_genes_tested - number_up_genes, 
                k = number_gene_in_gs, 
                lower.tail = FALSE)
  
  # getting UP genes within gene_set
  gene_list <- paste(eve_UP$gene.id[eve_UP$gene.id %in% gs$ensembl_gene], collapse = ",")
  # saving output
  gs_out <- rbind(gs_out, 
                 data.frame(gs_cat = unique(gs$gs_cat), gs_subcat = unique(gs$gs_subcat), gs_name = GS_name, gs_size = number_gene_in_gs, gs_updo = number_up_genes_in_gs, pval = Pval, shift = "UP", gene_names = gene_list))
  
  # Perform hypergeometric test DO
  number_do_genes <- nrow(eve_DO)
  number_do_genes_in_gs <- sum(eve_DO$gene.id %in% gs$ensembl_gene)
  
  Pval <- phyper(q = number_do_genes_in_gs - 1, 
                m = number_do_genes, 
                n = number_genes_tested - number_do_genes, 
                k = number_gene_in_gs, 
                lower.tail = FALSE)
  
  # getting DOWN genes within gene_set
  gene_list <- paste(eve_DO$gene.id[eve_DO$gene.id %in% gs$ensembl_gene], collapse = ",")
  # saving outout
  gs_out <- rbind(gs_out,
                  data.frame(gs_cat = unique(gs$gs_cat), gs_subcat = unique(gs$gs_subcat), gs_name = GS_name, gs_size = number_gene_in_gs, gs_updo = number_do_genes_in_gs, pval = Pval, shift = "DOWN", gene_names = gene_list))
}

gs_out <- gs_out %>% 
  mutate(pval_adj = p.adjust(pval, method = "BH")) %>%
  ungroup()

#write.table(gs_out, "output/GO_BP_Heart_Hgla.tab", quote = F, sep = "\t", row.names = F)
