# SCAVENGE

### Overview:

Co-localization approaches (such as gchromVAR) using genetic variants and single-cell epigenomic data are unfortunately uninformative for many cells given the extensive sparsity across single-cell profiles. Therefore, only a few cells from the truly relevant population demonstrate reliable phenotypic relevance. Nonetheless, the global high-dimensional features of individual single cells are sufficient to represent the underlying cell identities or states, which enables the relationships among such cells to be readily inferred15. By taking advantage of these attributes, SCAVENGE identifies the most phenotypically-enriched cells by co-localization and explores the transitive associations across the cell-to-cell network to assign each cell a probability representing the cell’s relevance to those phenotype-enriched cells via network propagation.

To address (2), we developed a novel enrichment method (**SCAVENGE**) (Single Cell Analysis of Variant Enrichment through Network propagation of GEnomic data) that can discriminate between closely related cell types and score single cells for GWAS enrichment. 




We've implemented **SCAVENGE** as an `R` package for computing single-cell based GWAS enrichments from fine-mapped posterior probabilities and quantitative epigenomic data (i.e. scATAC-seq and potentially other single-cell epigenome profiling).

### Installation:

Once all of the dependencies for `gchromVAR` are installed, the package can be installed 
directly from GitHub by typing the following into an `R` console:

```
devtools::install_github("https://github.com/sankaranlab/SCAVENGE")
```
### Tutorial:
This web resource and vignette compiliation shows how to reproduce these results in hematopoesis and how to run **SCAVENGE** on other data sets. 


### Citation:


### Contact:


