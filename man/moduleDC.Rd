% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moduleDC.R
\name{moduleDC}
\alias{moduleDC}
\title{Calculate modular differential connectivity (MDC)}
\usage{
moduleDC(inputMat = inputMat, inputMatB = NULL, design = design,
  compare = compare, genes = genes, labels = labels,
  gene_avg_signif = 0.05, number_DC_genes = 3, adjust = "bonferroni",
  bySlope = FALSE, corrType = "bicor", parallel = F,
  batchConfig = system.file("config/batchConfig_Local.R", package =
  "superNOVA"), split)
}
\arguments{
\item{inputMat}{The matrix (or data.frame) of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.}

\item{inputMatB}{The second matrix (or data.frame) to differentially correlate with the first of values (e.g., gene expression values from an RNA-seq or microarray study) that you are interested in analyzing. The rownames of this matrix should correspond to the identifiers whose correlations and differential correlations you are interested in analyzing, while the columns should correspond to the rows of the design matrix and should be separable into your groups.}

\item{design}{A standard model.matrix created design matrix. Rows correspond to samples and colnames refer to the names of the conditions that you are interested in analyzing. Only 0's or 1's are allowed in the design matrix. Please see vignettes for more information.}

\item{compare}{Vector of two character strings, each corresponding to one group name in the design matrix, that should be compared.}

\item{genes}{A character vector specifying gene symbols, present as rows in the inputMat, corresponding to each module label in the labels argument.}

\item{labels}{A character vector specifying module label names, one for each gene symbol in the genes argument, with overlap allowed (i.e., each gene can be in more than one module).}

\item{gene_avg_signif}{Minimum average gene level significance to report in the module.}

\item{number_DC_genes}{The number of differentially GOC or LOC genes to output in the final result.}

\item{adjust}{The method to adjust p-values (adjusts for multiple modules). Same inputs as p.adjust().}

\item{bySlope}{Whether or not to use the slopes instead of the correlations when calling GOC or LOC in the output for the two group case.}

\item{corrType}{The correlation type out of Pearson, Spearman, and bicor for the underlying linear model.}

\item{parallel}{Run the pipeline as single core or multithreaded? True or False.}

\item{batchConfig}{When running in parallel, the path to the batchtools configuration file. Best to run locally.}

\item{split}{When running in parallel, the number of times to split the input.}
}
\value{
A data frame with the module labels, the average change in difference in z-score between conditions (i.e., one measure of the modular average differential connectivity, or MeDC), and the empirical p-value for the significance of the change in correlation.
}
\description{
Takes modules of genes (possibly overlapping) and calculates the change in correlation among those genes between two conditions. Also reports the genes with the strongest gain in connectivity (i.e., average difference in z-score of > 0) and the strongest loss of correlation between conditions for each module, if any pass the significance measure specified.
}
\examples{
data(design_mat)
data(darmanis)
rownames(design_mat) = colnames(darmanis)
module_genes = list(mod1 = rownames(darmanis)[1:100],
 mod2 = rownames(darmanis)[90:190], mod3 = rownames(darmanis)[190:290])
modules = stack(module_genes)
modules$ind = as.character(modules$ind)
moduleDC_res = moduleDC(inputMat = darmanis, design = design_mat,
 compare = c("oligodendrocyte", "neuron"), genes = modules$values,
 labels = modules$ind)
}
