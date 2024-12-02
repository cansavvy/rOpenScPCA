This directory contains scripts to download and preprocess data used within the rOpenScPCA package.
There are also notebooks that explore some of the prepared datasets for exploration.

## Building Gene References

- The `build_gene_references.R` script creates a table of Ensembl id to gene symbol references.
It starts with the example ScPCA-formatted SCE object stored in test data, extracting the table of gene ids and gene symbols.
This is combined with the reference information extracted from example 10x Genomics datasets.
The full table is saved in `data/scpca_gene_reference.rda`.

- The `explore_gene_references.Rmd` notebook explores the resulting gene references table a bit to see where some of the conversions differ among different references.
