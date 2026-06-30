# README — RNA-seq counts, EVE/DESeq2 results, and orthologs (2025-09-12)

## Experimental design
This dataset includes bulk RNA-seq data from heart and liver tissues in four rodent species:  
- Naked mole-rat (Hgla, *Heterocephalus glaber*)  
- Damaraland mole-rat (Fdam, *Fukomys damarensis*)  
- Guinea pig (Cpor, *Cavia porcellus*)  
- Mouse (Mmus, *Mus musculus*)  

### Data sources
- **This study** — newly generated RNA-seq libraries.  
- **Berthelot et al. 2018** - https://doi.org/10.1038/s41559-017-0377-2  
- **Bens et al. 2018** - https://doi.org/10.1186/s12915-018-0546-z  
- **Faulkes et al. 2024** - https://doi.org/10.1038/s41467-024-46470-x  

---

## Contents of the folder
- `Heart_salmon_matrix_raw_counts.tab`
- `Heart_salmon_matrix_normalised_counts.tab`
- `Liver_salmon_matrix_raw_counts.tab`
- `Liver_salmon_matrix_normalised_counts.tab`
- `combined_deseq2Xeve_results.tab`
- `one2one_4species_EnsemblXZoonomia.tab`
- `species_tree.nwk`  

---

## How the tables were generated (short methods)

### RNA-seq processing & quantification
- **QC & trimming:** standard adapter/quality trimming.
- **Quantification:** Salmon pseudo-mapping at the transcript level against each species’ reference transcriptome.
- **Gene-level summarisation:** tximport to aggregate transcript estimates to genes with effective-length correction.  

### Normalisation for EVE
- **Normalised matrices** (`*_normalised_counts.tab`): counts were **scaled by effective transcript length** (as estimated by Salmon, txiimport+DESeq2), then normalised for library size (DESeq2 size factors) followed by a log2(counts + 1) transformation.  

### Phylogenetic expression shifts (EVE)
- **EVE** tests branch-specific expression shifts on the species tree (see branches bellow).
- For each gene × tissue × branch, EVE outputs a shift statistic and empirical p-value (permutation-based).

### Differential expression (DESeq2)
- Per tissue, DESeq2 run on tximport gene counts with size-factor & dispersion estimation.
- **Contrasts** align with the branch/species used for EVE (see bellow).
- Reported statistics in `combined_deseq2Xeve_results.tab` are from the Wald test (log2FC, SE, p-value, FDR).

### Branch design
- `Fdam` = **Damaraland mole-rat** (*Fukomys damarensis*)
- `Hgla` = **Naked mole-rat** (*Heterocephalus glaber*)
- `MR`   = **Ancestral mole-rat branch** (Naked + Damara mole-rat)

### Ortholog mapping
- `one2one_4species_EnsemblXZoonomia.tab` lists the ensembl ID of 1:1 orthologous genes across the four species.
- annotation source: Ensembl 102 + Zoonomia

---

## File summaries and column descriptions

### 1) `Heart_salmon_matrix_raw_counts.tab`
Gene-level raw counts (sum of Salmon `NumReads` across transcripts per gene via tximport) for heart samples.  
**Rows:** Mouse Ensembl 102 gene IDs.  
**Columns:** Samples, named `<Species>_<Tissue>_<SampleID>` (e.g., `Cpor_Heart_Rep1`, `Fdam_Heart_faulkes2`, `Hgla_Heart_Rep3`, `Mmus_Heart_faulkes5`).  

### 2) `Liver_salmon_matrix_raw_counts.tab`
Same format as above, for **liver**.

### 3) `Heart_salmon_matrix_normalised_counts.tab`
Gene-level normalised counts for heart samples.   
**Rows/Columns:** Same rows/columns as heart raw counts.  

### 4) `Liver_salmon_matrix_normalised_counts.tab`
Same format as above, for **liver**.

### 5) `combined_deseq2Xeve_results.tab`
Per gene × tissue × branch summary combining EVE branch-shift output with DESeq2 results aligned to the same branch labels.  
**Columns:**
- `gene.id` — Mouse Ensembl gene ID.
- `tissue` — `Heart` or `Liver`.
- `branch` — `Fdam`, `Hgla`, or `MR` (ancestral mole-rat).
- `GeneName` — Gene symbol (if available).
- **EVE:** `shift` (branch shift estimate), `empirical.pval`, `gene_shift` (discrete call: `UP`/`DOWN`/`N`).
- **DESeq2:** `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`, `gene_de` (direction call), `status` (`DE`/`N`).
- **Integration:** `gene_shift2` — qualitative concordance summary (e.g., `Both`, `DOWN`, `N`).

### 6) `one2one_4species_EnsemblXZoonomia.tab`
List of 1:1 ortholog genes included in the study with Ensembl gene ID concordance accross species.  
**Columns:** `GeneName`, `Mmus`, `Cpor`, `Fdam`, `Hgla` (Ensembl gene IDs).  

### 7) `species_tree.nwk`  
Newick species tree used for EVE modelling.  

---

## Contact
Questions on anything: Maëlle Daunesse (mdauness@pasteur.fr).