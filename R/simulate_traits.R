#' Simulate Quantitative Traits from PLINK Genotypes
#'
#' This function simulates a quantitative trait based on additive and epistatic 
#' genetic effects using genotype data from a PLINK dataset. The simulated trait 
#' is saved to a specified output file in a phenotype format compatible with
#' PLINK.
#'
#' @param plink_file Character. Path to the PLINK dataset (without file
#' extension). The function will append `.bed`, `.bim`, and `.fam` extensions
#' as needed.
#' @param output_file Character. Path to the output file where the simulated
#' trait will be saved.
#' @param additive_heritability Numeric. A value between 0 and 1 specifying the
#' proportion of trait variance due to additive genetic effects.
#' @param gxg_heritability Numeric. A value between 0 and 1 specifying the
#'  proportion of trait variance due to gene-by-gene (epistatic) interactions.
#'  The sum of `additive_heritability` and `gxg_heritability` must not exceed 1.
#' @param additive_indices Integer vector. Indices of SNPs contributing to
#' additive genetic effects.
#' @param gxg_indices_1 Integer vector. Indices of SNPs in the first group for
#' epistatic interactions.
#' @param gxg_indices_2 Integer vector. Indices of SNPs in the second group for
#' epistatic interactions.
#' @param log_level Character. Logging level for messages
#' (e.g., "DEBUG", "INFO", "WARNING"). Default is "WARNING".
#'
#' @return None. The simulated trait is written to the specified `output_file`.
#'
#' @details
#' The function uses the following components to simulate the trait:
#' - Additive genetic effects: Determined by `additive_indices` and the
#'   specified `additive_heritability`.
#' - Epistatic interactions: Simulated using pairs of SNPs from `gxg_indices_1`
#'   and `gxg_indices_2`, contributing to the `gxg_heritability`.
#' - Environmental effects: Any remaining variance not explained by genetic
#'   effects is assigned to random environmental noise.
#' 
#' The output file is in PLINK-compatible phenotype format with three columns: 
#' Family ID (`FID`), Individual ID (`IID`), and the simulated trait (`TRAIT`).
#'
#' @examples
#' plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package = "mmer"))
#' out_file <- tempfile()
#' additive_heritability <- 0.3
#' gxg_heritability <- 0.1
#' additive_snps <- sort(sample(1:100, 50, replace = FALSE))
#' gxg_group_1 <- sort(sample(additive_snps, 10, replace = FALSE))
#' gxg_group_2 <- sort(sample(setdiff(additive_snps, gxg_group_1), 10, replace = FALSE))
#' n_samples <- 200
#' simulate_traits(
#'   plink_file,
#'   out_file,
#'   additive_heritability,
#'   gxg_heritability,
#'   additive_snps,
#'   gxg_group_1,
#'   gxg_group_2
#' )
#' from_file <- read.table(out_file, header = TRUE)
#' head(from_file)
#'
#' @useDynLib mmer
#' @import genio
#' @import dplyr
#' @importFrom utils write.table
#' @export
simulate_traits <- function(plink_file,
                            output_file,
                            additive_heritability,
                            gxg_heritability,
                            additive_indices,
                            gxg_indices_1,
                            gxg_indices_2,
                            log_level = "WARNING") {


  logging::basicConfig(level = log_level)
  log <- logging::getLogger("mmer::simulate_traits")

  if (additive_heritability + gxg_heritability > 1) {
    stop("Additive heritability and gxg heritability should sum to less than 1")
  } else if (additive_heritability < 0 || gxg_heritability < 0) {
    stop("Heritabilities should be positive")
  }

  sim <- simulate_traits_cpp(
    plink_file,
    additive_heritability,
    gxg_heritability,
    additive_indices - 1,
    gxg_indices_1 - 1,
    gxg_indices_2 - 1
  )
  log$info(
    "Simulated additive variance %.2f, gxg variance %.2f, error variance %.2f",
    sim$additive_variance,
    sim$gxg_variance,
    sim$error_variance
  )
  fam_data <- read_fam(paste0(plink_file, ".fam"), verbose = FALSE)
  pheno_data <- data.frame(
    FID = fam_data$fam,
    IID = fam_data$id,
    TRAIT = sim$trait - mean(sim$trait)
  )
  write.table(
    pheno_data,
    file = output_file,
    sep = " ",
    quote = FALSE,
    row.names = FALSE
  )
}
