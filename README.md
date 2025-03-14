Overview

This repository contains code and analysis steps for an RNA-seq study designed to investigate gene expression dynamics in different cellular states. Specifically, the study compares gene expression between three distinct IMR90 fibroblast cell types:

Proliferating Cells (Proliferating)
Replicative Senescent Cells (Senescent)
Replicative Senescent Cells with Depleted Mitochondria (Mitochondria-Depleted Senescent)
The analysis involves various statistical and bioinformatics methods, including Principal Component Analysis (PCA), Volcano Plots, Heatmaps, Pathway Analysis, and Gene Expression Comparison, to better understand the molecular signatures of these different cell types and how mitochondrial depletion affects cellular senescence.


The purpose of this analysis is to:

Identify and compare gene expression patterns between proliferating cells, senescent cells, and mitochondria-depleted senescent cells.
Investigate the role of mitochondria depletion in modulating cellular senescence and the potential compensatory mechanisms.
Characterize molecular pathways implicated in senescence and mitochondrial dysfunction, focusing on changes in immune responses, metabolic adaptations, and cellular morphology.
Generate visualizations to help interpret complex gene expression data and discover underlying biological mechanisms.

Steps and Methodology

The analysis consists of several stages, each designed to answer a specific research question:

1. Data Quality Control and Preprocessing

Initial quality checks on the RNA-seq data were performed using MA plots, which helped in identifying the dispersion of gene expression across the different sample groups.
The analysis focused on normalizing the data and removing potential outliers to ensure the accuracy of downstream analyses.

2. Principal Component Analysis (PCA)

PCA was performed to identify clustering patterns and variance across the three cell types. This analysis revealed distinct groupings between proliferating, senescent, and mitochondria-depleted senescent cells, with the first two principal components explaining over 71% of the variance in the data.

3. Differential Gene Expression Analysis

Volcano plots were generated to identify the most significant genes upregulated or downregulated in the comparisons between:
Senescent vs. Proliferating (SP)
Senescent vs. Mitochondria-depleted Senescent (MTS)
Mitochondria-depleted Senescent vs. Proliferating (MTP)
A box plot of the top 15 most significant genes was produced to illustrate the differences in gene expression levels across these groups.

4. Heatmap and Clustering Analysis

Heatmaps were generated to visually display the clustering of significant genes across the different sample types. This allowed for the identification of gene expression patterns and trends unique to each experimental condition.

5. Pathway Analysis

Pathway analysis using Over-Representation Analysis (ORA) was conducted to identify the top pathways associated with upregulated and downregulated genes:
Upregulated genes in senescent cells were linked to immune response and extracellular matrix organization.
Downregulated genes in mitochondria-depleted senescent cells were related to extracellular matrix organization.
The results suggested significant shifts in metabolism and cellular response pathways, especially in the context of mitochondrial dysfunction.

6. Correlation and Signature Analysis

Pearson correlation coefficients were computed to assess the relationship between the different experimental groups. These showed positive correlations between MTP and SP, and MTP and MTS, while a negative correlation was found between MTS and SP.
Signature analyses, including cnetplots, helped identify the main genes and their associated pathways for each condition.
