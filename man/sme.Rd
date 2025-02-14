% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sme.R
\name{sme}
\alias{sme}
\title{Sparse Marginal Epistasis Test (SME)}
\usage{
sme(
  plink_file,
  pheno_file,
  mask_file = NULL,
  gxg_indices = NULL,
  chunk_size = NULL,
  n_randvecs = 10,
  n_blocks = 100,
  n_threads = 1,
  gxg_h5_group = "gxg",
  ld_h5_group = "ld",
  rand_seed = -1,
  log_level = "WARNING"
)
}
\arguments{
\item{plink_file}{Character. File path to the PLINK dataset
(without *.bed extension).
The function will append \code{.bim}, \code{.bed}, and \code{.fam} extensions
automatically.
The genotype data must not have any missing genotypes. Use PLINK to remove
variants with missing genotypes or impute them.}

\item{pheno_file}{Character. File path to a phenotype file in PLINK format.
The file should contain exactly one phenotype column.}

\item{mask_file}{Character or NULL. File path to an HDF5 file specifying
per-SNP masks for gene-by-gene interaction tests. This file informs which
SNPs are tested for marginal epistasis. Defaults to \code{NULL}, indicating no
masking. Masking impacts the scaling of memory and time.}

\item{gxg_indices}{Integer vector or NULL. List of indices corresponding to
SNPs to test for marginal epistasis.
If \code{NULL}, all SNPs in the dataset will be tested.
These indices are \strong{1-based}.}

\item{chunk_size}{Integer or NULL. Number of SNPs processed per chunk.
This influences memory
usage and can be left \code{NULL} to automatically determine the chunk size
based on \code{gxg_indices} and number of threads.}

\item{n_randvecs}{Integer. Number of random vectors used for stochastic trace
estimation.
Higher values yield more accurate estimates but increase computational
cost. Default is 10.}

\item{n_blocks}{Integer. Number of blocks into which SNPs are divided for
processing.
This parameter affects memory requirements. Default is 100.}

\item{n_threads}{Integer. Number of threads for OpenMP parallel processing.
Default is 1.}

\item{gxg_h5_group}{Character. Name of the HDF5 group within the mask file
containing gene-by-gene
interaction masks. SNPs in this group will be included in the gene-by-gene
interactions. Defaults to "gxg".}

\item{ld_h5_group}{Character. Name of the HDF5 group within the mask file
containing linkage disequilibrium
masks. SNPs in this group are excluded from analysis. Defaults to "ld".}

\item{rand_seed}{Integer. Seed for random vector generation. If \code{-1}, no seed
is set. Default is -1.}

\item{log_level}{Character. Logging level for messages. Must be in uppercase
(e.g., "DEBUG", "INFO", "WARNING", "ERROR"). Default is "WARNING".}
}
\value{
A list containing:
\itemize{
\item \code{summary}: A tibble summarizing results for each tested SNP, including:
\itemize{
\item \code{id}: Variant ID.
\item \code{index}: Index of the SNP in the dataset.
\item \code{chromosome}: Chromosome number.
\item \code{position}: Genomic position of the SNP.
\item \code{p}: P value for the gene-by-gene interaction test.
\item \code{pve}: Proportion of variance explained (PVE) by gene-by-gene interactions.
\item \code{vc}: Variance component estimate.
\item \code{se}: Standard error of the variance component.
}
\item \code{pve}: A long-format tibble of PVE for all variance components.
\item \code{vc_estimate}: A long-format tibble of variance component estimates.
\item \code{vc_se}: A long-format tibble of standard errors for variance components.
\item \code{average_duration}: Average computation time per SNP.
}
}
\description{
SME fits a linear mixed model in order to test for marginal epistasis. It
concentrates the scans for epistasis to regions of the genome that have known
functional enrichment for a trait of interest.
}
\details{
This function integrates PLINK-formatted genotype and phenotype data to
perform marginal epistasis tests on a set of SNPs. Using stochastic trace
estimation, the method computes variance components for gene-by-gene
interaction and genetic relatedness using the MQS estimator. The process is
parallelized using OpenMP when \code{n_threads > 1}.

The memory requirements and computation time scaling can be optimized through
the parameters \code{chunk_size}, \code{n_randvecs}, and \code{n_blocks}.

\strong{Mask Format Requirements}

The mask file format is an HDF5 file used for storing index data for
the masking process. This format supports data retrieval by index.
Below are the required groups and datasets within the HDF5 file:

The required group names can be configured as input parameters.
The defaults are described below.
\itemize{
\item \strong{Groups}:
\itemize{
\item \code{ld}: Stores SNPs in LD with the focal SNP. These SNPs will be \strong{excluded}.
\item \code{gxg}: Stores indices of SNPs that the marginal epistasis test is conditioned on. These SNPs will be \strong{included}.
}
\item \strong{Datasets}:
\itemize{
\item \verb{ld/<j>}: For each focal SNP \verb{<j>}, this dataset contains indices of SNPs
in the same LD block as that SNP. These SNPs will be \strong{excluded} from the gene-by-gene interaction covariance matrix.
\item \verb{gxg/<j>}: For each focal SNP \verb{<j>}, this dataset contains indices of SNPs to \strong{include} in the
the gene-by-gene interaction covariance matrix for focal SNP \verb{<j>}.
}
}

\strong{Important}: All indices in the mask file data are \strong{zero-based}, matching the zero-based indices of the PLINK \code{.bim} file.
}
\examples{
plink_file <- gsub("\\\\.bed", "", system.file("testdata", "test.bed", package="smer"))
pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="smer")
mask_file <- ""

# Parameter inputs
chunk_size <- 10
n_randvecs <- 10
n_blocks <- 10
n_threads <- 1

# 1-based Indices of SNPs to be analyzed
n_snps <- 100
snp_indices <- 1:n_snps

sme_result <- sme(
  plink_file,
  pheno_file,
  mask_file,
  snp_indices,
  chunk_size,
  n_randvecs,
  n_blocks,
  n_threads
)
head(sme_result$summary)

}
\references{
Stamp, J., Pattillo Smith, S., Weinreich, D., & Crawford, L. (2025).
Sparse modeling of interactions enables fast detection of genome-wide
epistasis in biobank-scale studies. bioRxiv, 2025.01.11.632557.

Stamp, J., DenAdel, A., Weinreich, D., & Crawford, L. (2023).
Leveraging the genetic correlation between traits improves the detection of
epistasis in genome-wide association studies.
G3: Genes, Genomes, Genetics, 13(8), jkad118.

Crawford, L., Zeng, P., Mukherjee, S., & Zhou, X. (2017).
Detecting epistasis with the marginal epistasis test in genetic mapping
studies of quantitative traits. PLoS genetics, 13(7), e1006869.
}
