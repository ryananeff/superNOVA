#' @title Gene ontology of differential correlation-classified genes.
#' @description Extracts a data frame of the top enriched gene sets in gene ontology databases using the hypergeometric test for gene synmols that are members of gene pairs in each of the classes specified in the differentially correlated gene pairs input table. Default parameter settings are to take in a result table with HGNC symbols and convert them to Ensembl symbols for gene ontology testing.
#' @param superNOVA_res The table of differential correlations outputted from superNOVA. Expected to have pValDiff or pValDiff_adj columns as well as zScoreDiff, Gene1, +/- Classes columns.
#' @param universe Character vector of gene symbols which should be used as the background in the hypergeomtric test. If using this in the context of a superNOVA experiment, this gene list most likely should be the gene set post-filtering, but prior to differential correlation analysis.
#' @param pval_gene_thresh p-value threshold to call a gene as having significant differential correlation or not.
#' @param classes Logical indicator specifying whether individual differential correlation gene classes should be extracted from the table or not. If not, compare MUST be set. Default = TRUE
#' @param compare Names of one or two groups to compare. If set to one group name, compare that group to the global model. If set to two group names, pairwise compare with each other. Default = NULL.
#' @param bySlope Whether to use the slope instead of the correlation rho when comparing one or two groups. This parameter has no effect when classes=TRUE.
#' @param geneNameCol Character vector specifying the name of the columns that are used to extract the gene symbols. Note that the default is c("Gene1", "Gene2"), but this only makes sense in the context of a full DGCA experiment. In the case of a splitSet, you may want to use "Gene1" to avoid counting the splitSet names in all of the categories.
#' @param pval_GO_cutoff Cutoff for the unadjusted p-values of gene ontology terms in the enrichment tests that should be displayed in the resulting table.
#' @param adjusted Logical indicating whether adjusted p-values from the differential correlation table (i.e., column "pValDiff_adj", when adjusted = TRUE) or unadjusted p-values (i.e., column "pValDiff", when adjusted = FALSE) should be used to subset the table into significant and non-significant portions.
#' @param HGNC_clean Logical indicating whether the input gene symbols should be switched to clean HGNC symbols using the checkGeneSymbols function from the R package HGNChelper. Only applies if HGNC symbols are inputted.
#' @param HGNC_switch Logical indicating whether or not the input gene symbols need to be switched from HGNC to Ensembl, the latter of which is required for GOstats enrichment test. Note that this is done by selecting the first Enembl symbol that maps to a particular HGNC symbol, which is not always unique. If you need more precision on the conversion, you should do this outside of the function and insert the Ensembl list to the function.
#' @param gene_ontology A string specifying the branch of GO that should be used for enrichment analysis. One of "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all". If "all" is chosen, then this function finds the enrichment for all of the terms and combines them into one table. Default = "all"
#' @param unique_genes Logical, if TRUE indicates that unique gene symbols within gene pairs from each category compared to the other groups should be chosen prior to GO enrichment analysis.
#' @param annotation The library indicating the GO annotation database from which the Go terms should be mapped to gene symbols. Default = "org.Hs.eg.db", which is the table for Homo sapiens. Other common choices include "org.Mm.eg.db", "org.Rn.eg.db". The corresponding annotation library needs to be installed.
#' @param calculateVariance Optionally, find the variance of the odds ratio for each enrichment test. In particular, this finds the standard error of the log odds ratio, which converges to a normal distribution much more quickly than the non-log OR.
#' @param conditional Logical specifying whether the GO analysis should be done conditionally to take into account the hierarchical structure of the GO database in making sense of the gene set enrichments.
#' @param superNOVA_find_significant Logical specifying whether this enrichment analysis should be performed on the result of a superNOVA analysis. If FALSE, then a superGO_res object, which is a named list of gene vectors, must be defined instead.
#' @param superGO_res Optional named list of gene vectors to find the enrichment of if superNOVA_find_signficiant is FALSE.
#' @import GOstats HGNChelper AnnotationDbi
#' @return A list of data frames corresponding to the gene ontology enrichment analysis results for the extracted gene sets from each of the differential correlation classes.
#' @references Agresti A: Categorical Data Analysis. 2012:70-77.
#' @export
superGO <- function(superNOVA_res, universe, pval_gene_thresh = 0.05,
  classes = TRUE, geneNameCol = c("Gene1", "Gene2"), pval_GO_cutoff = 1,
  HGNC_clean = TRUE, HGNC_switch = TRUE, gene_ontology = "all", adjusted = FALSE,
  annotation = "org.Hs.eg.db", conditional = FALSE, calculateVariance = FALSE,
  unique_genes = FALSE, regcor = FALSE, superNOVA_find_significant = TRUE, superGO_res = NULL,
  compare=NULL, bySlope=FALSE){

  if (!requireNamespace("GOstats", quietly = TRUE)) {
    stop("The R package GOstats is needed for the superGO function to work. Please install it.",
      call. = FALSE)
  }

  if(HGNC_clean){
    if (!requireNamespace("HGNChelper", quietly = TRUE)) {
      stop("The R package HGNChelper is needed for cleaning HGNC symbols with the HGNC_clean option you chose. Please install it.",
        call. = FALSE)
    }
  }

  if(HGNC_switch){
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop("The R package AnnotationDbi is needed for the superGO function to work with the annotation argument you chose. Please install it.",
        call. = FALSE)
    }
  }

  ##############################
  #set SAF to FALSE while restoring to default when the function is finished
  SAF = getOption("stringsAsFactors")
  on.exit(options(stringsAsFactors = SAF))
  options(stringsAsFactors = FALSE)

  if(!superNOVA_find_significant & is.null(superGO_res)){
    stop("If superNOVA_find_significant is FALSE, then superGO_res cannot be NULL; you must input a named list of gene vectors instead.")
  }

  if(HGNC_clean){
    universe = switchGenesToHGCN(universe)
  }
  if(HGNC_switch){
    univmap = AnnotationDbi::select(get(annotation), universe, "ENTREZID", "SYMBOL")
    universe = univmap[!duplicated(univmap[,1]), 2]
  }
  universe = universe[!is.na(universe)]
  universe = universe[!duplicated(universe)]

  if(superNOVA_find_significant){
    superGO_res = superNOVAFindSignificant(superNOVA_res, pval_gene_thresh = pval_gene_thresh,
      classes = classes, adjusted = adjusted, geneNameCol = geneNameCol,
      unique_genes = unique_genes, compare=compare,bySlope=bySlope)
    superGO_list = superGO_res
  } else {
    superGO_list = list()
  }
  for(i in 1:length(superGO_res)){

    gene_vector = superGO_res[[i]]

    if(length(gene_vector) == 0){
      name_df = paste0("enrichment_", names(superGO_res)[i])
      superGO_list[[name_df]] = list()
      next
    }

    if(HGNC_clean){
      gene_vector = switchGenesToHGCN(gene_vector)
    }
    if(HGNC_switch){
      univmap = AnnotationDbi::select(get(annotation), gene_vector, "ENTREZID", "SYMBOL")
      gene_vector = univmap[!duplicated(univmap[,1]), 2]
    }

    gene_vector = gene_vector[!is.na(gene_vector)]
    gene_vector = gene_vector[!duplicated(gene_vector)]

    superGO_df_tmp = findGOTermEnrichment(gene_vector = gene_vector,
      universe = universe,
      pval_GO_cutoff = pval_GO_cutoff, HGNC_clean = HGNC_clean,
      HGNC_switch = HGNC_switch, gene_ontology = gene_ontology,
      conditional = conditional, annotation = annotation)

    if(calculateVariance){
      if(length(gene_ontology) > 1 | gene_ontology == "all"){
        if(length(superGO_df_tmp) > 0){
          for(j in 1:length(superGO_df_tmp)){
            superGO_df_tmp_go = superGO_df_tmp[[j]]
            n11 = superGO_df_tmp_go$Count
            n21 = superGO_df_tmp_go$Size - superGO_df_tmp_go$Count
            n12 = superGO_df_tmp_go$gene_set_size - superGO_df_tmp_go$Count
            n22 = superGO_df_tmp_go$universe_size - n11 - n12 - n21
            superGO_df_tmp_go$LogOddsRatioSE = sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
            superGO_df_tmp[[j]] = superGO_df_tmp_go
          }
        }
      } else {
        n11 = superGO_df_tmp_go$Count
        n21 = superGO_df_tmp_go$Size - superGO_df_tmp_go$Count
        n12 = length(gene_vector) - superGO_df_tmp_go$Count
        n22 = length(universe) - n11 - n12 - n21
        superGO_df_tmp$LogOddsRatioSE = sqrt((1/n11) + (1/n12) + (1/n21) + (1/n22))
      }
    }
    name_df = paste0("enrichment_", names(superGO_res)[i])
    superGO_list[[name_df]] = superGO_df_tmp
  }
  return(superGO_list)

}

#' @title Switches a gene vector to cleaned HGNC symbols.
#' @description Where possible, switches a character vector of gene names to cleaned and updated HGNC symbols.
#' @param gene_list Character vector of gene names.
#' @return Character vector of cleaned gene names.
#' @export
switchGenesToHGCN <- function(gene_list){
  gene_listHGNC = suppressWarnings(HGNChelper::checkGeneSymbols(gene_list))
  to_change = which((gene_listHGNC$Approved == FALSE) &
    !is.na(gene_listHGNC$Suggested.Symbol))
  gene_list[to_change] = gene_listHGNC[to_change, ]$Suggested.Symbol
  gene_list = toupper(gene_list)
  return(gene_list)
}

#' @title Find GO enrichment for a gene vector (using GOstats).
#' @description Given a gene character vector and a universe character vector, which can be either Ensembl or HGNC symbols, find the over-representation enrichment of the gene list relative to the universe in a gene ontology category using the hypergeometric test and the GOstats R package.
#' @param gene_vector Character vector gene symbols of interest.
#' @param universe Character vector of gene symbols which should be used as the background in the hypergeomtric test. If using this in the context of a superNOVA experiment, this gene list most likely should be the gene set post-filtering, but prior to differential correlation analysis.
#' @param pval_GO_cutoff Cutoff for the p-values of gene ontology terms in the enrichment tests that should be displayed in the resulting table.
#' @param cleanNames Logical, indicating whether the gene names for the universe and gene vector should be cleaned prior to enrichment analysis.
#' @param HGNC_clean Logical indicating whether the input gene symbols should be switched to clean HGNC symbols using the checkGeneSymbols function from the R package HGNChelper. Only applies if HGNC symbols are inputted and cleanNames is TRUE.
#' @param HGNC_switch Logical indicating whether or not the input gene symbols need to be switched from HGNC to Ensembl, the latter of which is required for GOstats enrichment test. Note that this is done by selecting the first Enembl symbol that maps to a particular HGNC symbol, which is not always unique. If you need more precision, you should do this outside of the function and insert the Ensembl list to the function. Only applies if cleanNames is TRUE.
#' @param gene_ontology A string specifying the branch of GO that should be used for enrichment analysis. One of "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component), or "all". If "all" is chosen, then this function finds the enrichment for all of the terms and combines them into one table. Default = "all"
#' @param annotation The library indicating the GO annotation database from which the Go terms should be mapped to gene symbols. Default = "org.Hs.eg.db", which is the table for Homo sapiens. Other common choices include "org.Mm.eg.db", "org.Rn.eg.db". The corresponding annotation library needs to be installed.
#' @param conditional Logical specifying whether the GO analysis should be done conditionally to take into account the hierarchical structure of the GO database in making sense of the gene set enrichments. Default = TRUE.
#' @return A data frame with the term enrichments of the GO enrichment analysis given the input gene set and universe.
#' @export
findGOTermEnrichment <- function(gene_vector, universe, pval_GO_cutoff = 1, HGNC_switch = TRUE,
  HGNC_clean = TRUE, gene_ontology = "all", conditional = TRUE, annotation = "org.Hs.eg.db",
  cleanNames = FALSE){

  if(length(gene_vector) == 0){
    empty_list = list()
    return(empty_list)
  }

  if(cleanNames){
    if(HGNC_clean){
      gene_vector = switchGenesToHGCN(gene_vector)
    }
    if(HGNC_switch){
      genemap = AnnotationDbi::select(get(annotation), gene_vector, "ENTREZID", "SYMBOL")
      gene_vector = genemap[!duplicated(genemap[,1]), 2]
    }
    gene_vector = gene_vector[!is.na(gene_vector)]
    gene_vector = gene_vector[!duplicated(gene_vector)]
  }

  if(!gene_ontology == "all"){
    params = new("GOHyperGParams",
      geneIds = gene_vector,
      universeGeneIds = universe,
      annotation = annotation,
      ontology = gene_ontology,
      pvalueCutoff = pval_GO_cutoff,
      conditional = conditional,
      testDirection = "over")
    hgOver = GOstats::hyperGTest(params)
    hg_df = GOstats::summary(hgOver)
  }

  if(gene_ontology == "all"){
    params_bp = new("GOHyperGParams",
      geneIds = gene_vector,
      universeGeneIds = universe,
      annotation = annotation,
      ontology = "BP",
      pvalueCutoff = pval_GO_cutoff,
      conditional = conditional,
      testDirection = "over")
    hgOver_bp = GOstats::hyperGTest(params_bp)
    hg_df_bp = GOstats::summary(hgOver_bp)
    hg_df_bp$Ontology = rep("BP", nrow(hg_df_bp))

    set_sizes_bp = capture.output(print(hgOver_bp))
    gene_set_bp = set_sizes_bp[grep("gene set size", set_sizes_bp)]
    gene_set_size_bp = as.numeric(trimws(strsplit(gene_set_bp, ":")[[1]][2]))
    universe_set_bp = set_sizes_bp[grep("universe size", set_sizes_bp)]
    universe_size_bp = as.numeric(trimws(strsplit(universe_set_bp, ":")[[1]][2]))
    hg_df_bp$gene_set_size = rep(gene_set_size_bp, nrow(hg_df_bp))
    hg_df_bp$universe_size = rep(universe_size_bp, nrow(hg_df_bp))

    params_mf = new("GOHyperGParams",
      geneIds = gene_vector,
      universeGeneIds = universe,
      annotation = annotation,
      ontology = "MF",
      pvalueCutoff = pval_GO_cutoff,
      conditional = conditional,
      testDirection = "over")
    hgOver_mf = GOstats::hyperGTest(params_mf)
    hg_df_mf = GOstats::summary(hgOver_mf)
    hg_df_mf$Ontology = rep("MF", nrow(hg_df_mf))

    set_sizes_mf = capture.output(print(hgOver_mf))
    gene_set_mf = set_sizes_mf[grep("gene set size", set_sizes_mf)]
    gene_set_size_mf = as.numeric(trimws(strsplit(gene_set_mf, ":")[[1]][2]))
    universe_set_mf = set_sizes_mf[grep("universe size", set_sizes_mf)]
    universe_size_mf = as.numeric(trimws(strsplit(universe_set_mf, ":")[[1]][2]))
    hg_df_mf$gene_set_size = rep(gene_set_size_mf, nrow(hg_df_mf))
    hg_df_mf$universe_size = rep(universe_size_mf, nrow(hg_df_mf))

    params_cc = new("GOHyperGParams",
      geneIds = gene_vector,
      universeGeneIds = universe,
      annotation = annotation,
      ontology = "CC",
      pvalueCutoff = pval_GO_cutoff,
      conditional = conditional,
      testDirection = "over")
    hgOver_cc = GOstats::hyperGTest(params_cc)
    hg_df_cc = GOstats::summary(hgOver_cc)
    hg_df_cc$Ontology = rep("CC", nrow(hg_df_cc))

    set_sizes_cc = capture.output(print(hgOver_cc))
    gene_set_cc = set_sizes_cc[grep("gene set size", set_sizes_cc)]
    gene_set_size_cc = as.numeric(trimws(strsplit(gene_set_cc, ":")[[1]][2]))
    universe_set_cc = set_sizes_cc[grep("universe size", set_sizes_cc)]
    universe_size_cc = as.numeric(trimws(strsplit(universe_set_cc, ":")[[1]][2]))
    hg_df_cc$gene_set_size = rep(gene_set_size_cc, nrow(hg_df_cc))
    hg_df_cc$universe_size = rep(universe_size_cc, nrow(hg_df_cc))

    hg_df = list(BP = hg_df_bp, MF = hg_df_mf, CC = hg_df_cc)
  }
  return(hg_df)
}
