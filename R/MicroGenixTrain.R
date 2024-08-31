
#' Training Models with Genotype-Microbiome Interactions
#'
#' @description MicroGenixTrain is used to fit models to predict gene expression with elastic net. Predictors include SNPs, microbes, and SNP-microbe interactions.
#'
#' @param input_taxa the input file name of taxonomic abundance profiles.
#' @param input_geno the input file name of host genotype dosage (0/1/2).
#' @param input_expr the input file name of host gene expression levels.
#' @param log_trans_taxa logical. If TRUE the taxonomic profiles are log2-transformed.
#' @param log_trans_taxa_add a numeric value of pseudo-count added to the taxonomic profiles to avoid taking the logarithm of zero.
#' @param log_trans_gene logical. If TRUE the gene expression levels are log2-transformed.
#' @param log_trans_gene_add a numeric value of pseudo-count added to the gene expression levels to avoid taking the logarithm of zero.
#' @param fold_id "foldid" for cv.glmnet, see \code{\link[glmnet]{cv.glmnet}}
#' @param num_folds "nfolds" for cv.glmnet, see \code{\link[glmnet]{cv.glmnet}}
#' @param output_prefix the prefix of output file names.
#' @param pval_fit logical. This is an experimental argument. If TRUE fit a basic linear model using variables selected by elastic net to obtain p-values.
#'
#' @return A model object and coefficients.
#' @export
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet coef.glmnet
#' @importFrom dplyr left_join
#' @importFrom stats as.formula
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats predict
#' @importFrom utils read.csv
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' input_taxa <- 'extdata/example_taxon_abundance.csv'
#' input_geno <- 'extdata/example_genotype_dosage.csv'
#' input_expr <- 'extdata/example_gene_expression.csv'
#' fit_data <- MicroGenixTrain(input_taxa, input_geno, input_expr)
#' }
MicroGenixTrain <- function(input_taxa, input_geno, input_expr,
                            log_trans_taxa = TRUE, log_trans_taxa_add = 0.01,
                            log_trans_gene = FALSE, log_trans_gene_add = 0.01,
                            fold_id = NULL, num_folds = 10, output_prefix = NULL,
                            pval_fit = FALSE){
  taxa <- utils::read.csv(input_taxa, check.names = F, row.names = 1)
  if(log_trans_taxa){
    taxa <- log2(taxa + log_trans_taxa_add)
  }
  geno <- utils::read.csv(input_geno, check.names = F, row.names = 1)
  expr <- utils::read.csv(input_expr, check.names = F, row.names = 1)
  if(log_trans_gene){
    expr <- log2(expr + log_trans_gene_add)
  }
  sample_names <- intersect(rownames(taxa), rownames(geno))
  sample_names <- intersect(sample_names, rownames(expr))
  num_samples <- length(sample_names)
  if(num_samples == 0){
    stop(paste0("Sample names in input files do not match"))
  }
  print(paste('Training on', num_samples, 'samples present in all input files'))
  taxa <- taxa[sample_names, , drop = FALSE]
  geno <- geno[sample_names, , drop = FALSE]
  expr <- expr[sample_names, , drop = FALSE]
  colnames(taxa) <- gsub(' ', '-', colnames(taxa), fixed = T)
  colnames(geno) <- gsub(' ', '-', colnames(geno), fixed = T)
  colnames(taxa) <- paste0('MICROBE_', colnames(taxa))
  colnames(geno) <- paste0('SNP_', colnames(geno))
  x <- as.matrix(cbind(taxa, geno))
  if(ncol(expr) > 1){
    stop("Modeling should be limited to one gene at a time")
  }
  y <- expr[,1]
  fit_no_interaction <- glmnet::cv.glmnet(x, y, alpha = 0.5, nfolds = num_folds, foldid = fold_id)
  df_coef <- as.data.frame(as.matrix(glmnet::coef.glmnet(fit_no_interaction, s="lambda.min")))
  var_nonzero <- rownames(df_coef)[df_coef[,1] != 0]
  var_taxa <- var_nonzero[grep('MICROBE_', var_nonzero, fixed = T)]
  var_geno <- var_nonzero[grep('SNP_', var_nonzero, fixed = T)]
  print(paste('The number of microbial variables with non-zero coefficients is', length(var_taxa)))
  print(paste('The number of SNP variables with non-zero coefficients is', length(var_geno)))
  if(length(var_taxa) > 0 & length(var_geno) > 0){
    pick_taxa <- taxa[, var_taxa, drop = FALSE]
    pick_geno <- geno[, var_geno, drop = FALSE]
    stopifnot(all(rownames(pick_taxa) == rownames(pick_geno)))
    mtx_interact <- matrix(nrow = nrow(pick_geno), ncol = length(var_geno)*length(var_taxa))
    rownames(mtx_interact) <- rownames(pick_geno)
    feature_interaction <- c()
    k <- 0
    for(i in 1:length(var_taxa)){
      for(j in 1:length(var_geno)){
        feature_interaction <- c(feature_interaction, paste0(var_taxa[i], '__x__', var_geno[j]))
        k <- k + 1
        mtx_interact[,k] <- pick_taxa[,i] * pick_geno[,j]
      }
    }
    colnames(mtx_interact) <- feature_interaction
    x <- cbind(mtx_interact, x)
    fit <- glmnet::cv.glmnet(x, y, alpha = 0.5, nfolds = num_folds, foldid = fold_id)
    if(!is.null(output_prefix)){
      saveRDS(fit, file = paste0(output_prefix, '.rds'))
    }
    df_coef <- as.matrix(glmnet::coef.glmnet(fit, s="lambda.min"))
    df_coef <- data.frame(feature = rownames(df_coef), coefficient = df_coef[,1])
    df_coef <- df_coef[df_coef[,1] != '(Intercept)' & df_coef[,2] != 0,]
    rownames(df_coef) <- NULL
    if(pval_fit){
      x_lm <- x[,colnames(x) %in% df_coef$feature]
      train_data <- as.data.frame(x_lm)
      train_data$express <- y
      fit_lm <- stats::lm(express ~ ., data = train_data)
      fit_lm_coef <- summary(fit_lm)$coefficients
      fit_lm_coef <- data.frame(feature = gsub('`','',rownames(fit_lm_coef)),
                                pvalue_lm = fit_lm_coef[, 'Pr(>|t|)'])
      df_coef <- dplyr::left_join(df_coef, fit_lm_coef, by = 'feature')
    }
    if(!is.null(output_prefix)){
      utils::write.table(df_coef, file = paste0(output_prefix, '.tsv'), row.names = F, quote = F, sep = '\t')
    }
    result <- list(model = fit, coefficient = df_coef)
  } else {
    print('Modeling has been omitted due to a lack of sufficient variables')
    result <- NULL
  }
  return(result)
}
