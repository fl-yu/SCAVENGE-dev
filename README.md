# SCAVENGE

### About:

Two outstanding challenges in the post-GWAS era are (1) the precise identification of causal variants within associated loci and (2) determining the exact mechanisms by which these variants result in the observed phenotypes, starting with identification of the pertinent cell types. To address (1), we used robust genetic fine mapping to identify hundreds of likely causal variants for 16 blood cell traits, allowing for up to 5 causal variants in each locus. We combined our fine-mapped results with high resolution open chromatin data for 18 primary hematopoietic populations and derived functional annotations to identify predicted target genes, mechanisms, and disease relevance. Moreover, we elucidate compelling anecdotes for the utility of this approach. To address (2), we developed a novel enrichment method (**gchromVAR**) that can discriminate between closely related cell types and score single cells for GWAS enrichment. 

We've implemented **gchromVAR** as an `R` package for computing cell-type specific GWAS enrichments from GWAS summary statistics and quantitative epigenomic data. This web resource and vignette compiliation shows how to reproduce these results in hematopoesis and how to run **gchromVAR** on other data sets. 

### Installation:

Once all of the dependencies for `gchromVAR` are installed, the package can be installed 
directly from GitHub by typing the following into an `R` console:

```
devtools::install_github("caleblareau/gchromVAR")
```
