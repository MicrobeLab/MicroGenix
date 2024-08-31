
#' Predicting Host Gene Expression Using Models with Genotype-Microbiome Interactions
#'
#' @description MicroGenixPredict predicts gene expression using models trained with MicroGenixTrain.
#'
#' @param model the file name of a model object (.rds).
#' @param input_taxa the input file name of taxonomic abundance profiles.
#' @param input_geno the input file name of host genotype dosage (0/1/2).
#' @param observed_expr the file name of observed (actual) gene expression levels. If provided, accuracy is evaluated.
#' @param log_trans_taxa logical. If TRUE the taxonomic profiles are log2-transformed.
#' @param log_trans_taxa_add a numeric value of pseudo-count added to the taxonomic profiles to avoid taking the logarithm of zero.
#' @param log_trans_gene logical. If TRUE the observed gene expression levels are log2-transformed.
#' @param log_trans_gene_add a numeric value of pseudo-count added to the observed gene expression levels to avoid taking the logarithm of zero.
#'
#' @return Predicted gene expression levels and R2 (if observed_expr is provided).
#' @export
#'
#' @importFrom glmnet coef.glmnet
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
#' input_taxa <- system.file('extdata', 'example_taxon_abundance.csv', package="MicroGenix")
#' input_geno <- system.file('extdata', 'example_genotype_dosage.csv', package="MicroGenix")
#' model <- 'example_output_model.rds'
#' predicted_data <- MicroGenixPredict(model, input_taxa, input_geno)
#' }
MicroGenixPredict <- function(model, input_taxa, input_geno, observed_expr = NULL,
                            log_trans_taxa = TRUE, log_trans_taxa_add = 0.01,
                            log_trans_gene = FALSE, log_trans_gene_add = 0.01){
  taxa <- utils::read.csv(input_taxa, check.names = F, row.names = 1)
  if(log_trans_taxa){
    taxa <- log2(taxa + log_trans_taxa_add)
  }
  geno <- utils::read.csv(input_geno, check.names = F, row.names = 1)

  sample_names <- intersect(rownames(taxa), rownames(geno))
  num_samples <- length(sample_names)
  if(num_samples == 0){
    stop("Sample names in input files do not match")
  }
  taxa <- taxa[sample_names, , drop = FALSE]
  geno <- geno[sample_names, , drop = FALSE]
  colnames(taxa) <- gsub(' ', '-', colnames(taxa), fixed = T)
  colnames(geno) <- gsub(' ', '-', colnames(geno), fixed = T)
  colnames(taxa) <- paste0('MICROBE_', colnames(taxa))
  colnames(geno) <- paste0('SNP_', colnames(geno))
  x <- as.matrix(cbind(taxa, geno))
  fit <- readRDS(model)
  df_coef <- as.matrix(glmnet::coef.glmnet(fit, s="lambda.min"))
  vars <- rownames(df_coef)
  var_interaction <- vars[grep('__x__', vars, fixed = T)]
  if(length(var_interaction) > 0){
    mtx_interact <- matrix(nrow = nrow(x), ncol = length(var_interaction))
    colnames(mtx_interact) <- var_interaction
    var_interaction <- lapply(var_interaction, function(x) strsplit(x, "__x__")[[1]])
    for(k in 1:length(var_interaction)){
      mtx_interact[,k] <- taxa[,var_interaction[[k]][1]] * geno[,var_interaction[[k]][2]]
      k <- k + 1
    }
    x <- cbind(mtx_interact, x)
  } else {
    print('Note: The model does not include an interaction term')
  }
  predicted_expr <- predict(fit, newx = x, s="lambda.min")
  if(!is.null(observed_expr)){
    observed_expr <- utils::read.csv(observed_expr, check.names = F, row.names = 1)
    if(log_trans_gene){
      observed_expr <- log2(observed_expr + log_trans_gene_add)
    }
    sample_names <- intersect(rownames(predicted_expr), rownames(observed_expr))
    observed_expr <- observed_expr[sample_names, , drop = FALSE]
    predicted_expr <- predicted_expr[sample_names, , drop = FALSE]
    r2 <- 1 - sum((predicted_expr[,1] - observed_expr[,1]) ^ 2) / sum((observed_expr[,1] - mean(observed_expr[,1])) ^ 2)
    result <- list(predicted_expr = predicted_expr[,1], observed_expr = observed_expr[,1], r2 = r2)
  } else {
    result <- list(predicted_expr = predicted_expr[,1], observed_expr = NULL, r2 = NULL)
  }
  return(result)
}
