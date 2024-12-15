#' @title Simulated Dataset for Genome-Wide Interaction Analysis
#' @description
#' `getting_started` is a simulated dataset created to demonstrate the use of 
#' the `mme()` function for genome-wide interaction analyses. It contains 
#' results from a simulated analysis involving additive genetic effects and 
#' gene-by-gene (GxG) interactions.
#'
#' @details
#' The dataset was generated as follows:
#' 
#' - **Genotype Simulation**:  
#'   Genotype data for 5000 individuals and 6,000 SNPs was simulated with
#'   synthetic allele counts.
#' 
#' - **Phenotype Simulation**:  
#'   Phenotypic values were simulated with an additive heritability of 0.3 and a
#'   GxG interaction heritability of 0.25. A set of 100 SNPs were selected for
#'   additive effects, and two groups of 5 SNPs each were used for GxG
#'   interactions.
#' 
#' - **PLINK-Compatible Files**:  
#'   The simulated data was saved in PLINK-compatible `.bed`, `.fam`,
#'   and `.bim` files.
#' 
#' - **Interaction Analysis**:  
#'   The `mme()` function was used to perform genome-wide interaction analyses
#'   on a subset of SNP indices, including the GxG SNP groups and 100 additional
#'   additive SNPs. Memory-efficient computation parameters 
#'   (e.g., `chun_ksize`, `n_randvecs`, and `n_blocks`) were applied.
#'
#' @format
#' A list with results from `mme()`, including the following components:
#' \describe{
#'   \item{`summary`}{A data frame summarizing the analysis results, including
#'   p-values for SNP associations (`p`).}
#'   \item{`pve`}{A data frame containing the per SNP variance component
#'   estimates normalized to phenotypic variance explained (PVE).}
#'   \item{`vc`}{A data frame containing the per SNP variance component
#'   estimates.}
#'   \item{`gxg_snps`}{A vector containing the indices of the SNPs assigned to
#'   have epistatic interactions in the trait simulations.}
#' }
#'
#' @section Key Parameters:
#' - **Additive Heritability**: 0.3  
#' - **GxG Heritability**: 0.25  
#' - **Number of Samples**: 5000  
#' - **Number of SNPs**: 6,000  
#' - **Selected Additive SNPs**: 100  
#' - **Selected GxG SNP Groups**:  
#'   - Group 1: 5 SNPs  
#'   - Group 2: 5 SNPs  
#'
#' @usage 
#' data("getting_started")
#'
#' @examples
#' data("getting_started")
#' head(getting_started$summary)
#'
#' @seealso
#' \link[mme]{mme}
#'
#' @keywords datasets
#' @source data-raw/getting_started.R
#' @import mvMAPIT
"getting_started"