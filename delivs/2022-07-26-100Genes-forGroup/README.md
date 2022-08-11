# README.md

## Preliminaries

Data in this repository is from subset analysis:
	- subset 1 (42 LOAD; 42 controls)
	- subset 2 (42 LOAD; 42 *different* controls)
	- LFDR: local false discovery rate. We use $\alpha=0.05$ as significance cutoff
	- Top 100 genes determined by...
		1. Taking all genes with at least 50 CpGs in [Transcription Start Site, Transcription Stop Site]. (Note that this counts CpGs in introns and exons.)
		2. Ordering by average (between two subsets) percent of CpGs that are significant. This is meant to correct for bias induced by gene size (i.e. large genes are expected to have more significanct CpGs by virtue of having more CpGs).
		3. Keeping top 100.

## File Descriptions

- DMPs-by-gene-comparison-top100.csv contains information on the number of differentially methylated positions (DMPs) after performing analysis on two subsets separatelty.
- comp-gene.png summarizes the DMPs
	- As usual, one CpG is represented by one point in the Chicago plot
	- Above the line at about 1.2 is significant (in subset 1) and hypermethylated
	- Below the line is significant (in subset 1).
	- A handful of genes in the list do not have plots due to errors in package dependencies
	- Points are colored by whether they achieve significance in subset 2 analysis. So in a "perfect" world, the orange points == points above the line, and blue points are the same as the points below the line.
- 2022-07-26-nDMPs-vs-nCpGs.html is an interactive plot that shows you where each gene lies in `Number of CpGs in gene` $\times$ `(average between two subsets) number of significance CpGs per gene`. Hover over points to see what the gene is called.