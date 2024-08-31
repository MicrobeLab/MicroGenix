# MicroGenix

MicroGenix integrates multi-omics data to identify the genotype-microbiome interactions in shaping molecular phenotype.

## Installation

The MicroGenix R package is easily installed from the GitHub repository:

    install.packages("devtools") 
    devtools::install_github("MicrobeLab/MicroGenix")

## Manual

The MicroGenix R package includes three main functions. Detailed documentation is available through the R help interface (`?MicroGenixTrain`, `?MicroGenixPredict`, `?MicroGenixAssociation`). Demo input files for running examples below are available in the `extdata` folder.

### MicroGenixTrain

`MicroGenixTrain` is used to fit models to predict gene expression with elastic net. Predictors include SNPs, microbes, and SNP-microbe interactions.

    input_taxa <- 'extdata/example_taxon_abundance.csv'
    input_geno <- 'extdata/example_genotype_dosage.csv'
    input_expr <- 'extdata/example_gene_expression.csv'
    fit_data <- MicroGenixTrain(input_taxa, input_geno, input_expr, output_prefix = 'example_output_model')

`input_taxa` and `input_geno` are comma-separated table files in shape [number_of_samples, number_of_taxa] and [number_of_samples, number_of_SNPs], respectively. `input_expr` is a file in shape [number_of_samples,] with the expression levels of a single gene.

### MicroGenixPredict

`MicroGenixPredict` predicts gene expression using models trained with `MicroGenixTrain`.

    input_taxa <- 'extdata/example_taxon_abundance.csv'
    input_geno <- 'extdata/example_genotype_dosage.csv'
    model <- 'example_output_model.rds'
    predicted_data <- MicroGenixPredict(model, input_taxa, input_geno)

`input_taxa` and `input_geno` are in the same format as the input files for `MicroGenixTrain`. `model` is a model file generated using `MicroGenixTrain`.

### MicroGenixAssociation

`MicroGenixAssociation` performs association tests between phenotype and gene expression levels predicted using `MicroGenixPredict`.

    pred_expr <- predicted_data$predicted_expr
    metadata <- 'extdata/example_metadata.csv'
    assoc_result <- MicroGenixAssociation(pred_expr, metadata, pheno = 'pheno')

`metadata` is a comma-separated file with phenotype and covariates.

## Other Resources

Bugs and difficulties in using MicroGenix are welcome on [the issue tracker](https://github.com/MicrobeLab/MicroGenix/issues).

The MicroGenix method is described in "Multi-omics integration unravels genotype-microbiome interactions shaping the conjunctival transcriptome". Models trained using the conjunctival multi-omics data are available [here](https://drive.google.com/file/d/1kS_3GkNiEIoaGfC72E0o7mQTUCHsmN_p/view?usp=drive_link).
