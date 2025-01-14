#' Sparse Marginal Epistasis Test (SME)
#'
#' SME fits a linear mixed model in order to test for marginal epistasis. It
#' concentrates the scans for epistasis to regions of the genome that have known
#' functional enrichment for a trait of interest.
#'
#' @references Stamp, J., Pattillo Smith, S., Weinreich, D., & Crawford, L. (2025).
#' Sparse modeling of interactions enables fast detection of genome-wide
#' epistasis in biobank-scale studies. bioRxiv, 2025.01.11.632557.
#' @references Stamp, J., DenAdel, A., Weinreich, D., & Crawford, L. (2023).
#' Leveraging the genetic correlation between traits improves the detection of
#' epistasis in genome-wide association studies.
#' G3: Genes, Genomes, Genetics, 13(8), jkad118.
#' @references Crawford, L., Zeng, P., Mukherjee, S., & Zhou, X. (2017).
#' Detecting epistasis with the marginal epistasis test in genetic mapping
#' studies of quantitative traits. PLoS genetics, 13(7), e1006869.
#'
#'
#' @param plink_file Character. File path to the PLINK dataset
#'   (without *.bed extension).
#'   The function will append `.bim`, `.bed`, and `.fam` extensions
#'   automatically.
#'   The genotype data must not have any missing genotypes. Use PLINK to remove
#'   variants with missing genotypes or impute them.
#' @param pheno_file Character. File path to a phenotype file in PLINK format.
#'   The file should contain exactly one phenotype column.
#' @param mask_file Character or NULL. File path to an HDF5 file specifying
#'   per-SNP masks for gene-by-gene interaction tests. This file informs which
#'   SNPs are tested for marginal epistasis. Defaults to `NULL`, indicating no
#'   masking. Masking impacts the scaling of memory and time.
#' @param gxg_indices Integer vector or NULL. List of indices corresponding to
#'   SNPs to test for marginal epistasis.
#'   If `NULL`, all SNPs in the dataset will be tested.
#'   These indices are **1-based**.
#' @param chunk_size Integer or NULL. Number of SNPs processed per chunk.
#'   This influences memory
#'   usage and can be left `NULL` to automatically determine the chunk size
#'   based on `gxg_indices` and number of threads.
#' @param n_randvecs Integer. Number of random vectors used for stochastic trace
#'   estimation.
#'   Higher values yield more accurate estimates but increase computational
#'   cost. Default is 10.
#' @param n_blocks Integer. Number of blocks into which SNPs are divided for
#'   processing.
#'   This parameter affects memory requirements. Default is 100.
#' @param n_threads Integer. Number of threads for OpenMP parallel processing.
#'   Default is 1.
#' @param gxg_h5_group Character. Name of the HDF5 group within the mask file
#'   containing gene-by-gene
#'   interaction masks. SNPs in this group will be included in the gene-by-gene
#'   interactions. Defaults to "gxg".
#' @param ld_h5_group Character. Name of the HDF5 group within the mask file
#'   containing linkage disequilibrium
#'   masks. SNPs in this group are excluded from analysis. Defaults to "ld".
#' @param rand_seed Integer. Seed for random vector generation. If `-1`, no seed
#'   is set. Default is -1.
#' @param log_level Character. Logging level for messages. Must be in uppercase
#'   (e.g., "DEBUG", "INFO", "WARNING", "ERROR"). Default is "WARNING".
#'
#' @return A list containing:
#'   - `summary`: A tibble summarizing results for each tested SNP, including:
#'       - `id`: Variant ID.
#'       - `index`: Index of the SNP in the dataset.
#'       - `chromosome`: Chromosome number.
#'       - `position`: Genomic position of the SNP.
#'       - `p`: P value for the gene-by-gene interaction test.
#'       - `pve`: Proportion of variance explained (PVE) by gene-by-gene interactions.
#'       - `vc`: Variance component estimate.
#'       - `se`: Standard error of the variance component.
#'   - `pve`: A long-format tibble of PVE for all variance components.
#'   - `vc_estimate`: A long-format tibble of variance component estimates.
#'   - `vc_se`: A long-format tibble of standard errors for variance components.
#'   - `average_duration`: Average computation time per SNP.
#'
#' @details
#' This function integrates PLINK-formatted genotype and phenotype data to
#' perform marginal epistasis tests on a set of SNPs. Using stochastic trace
#' estimation, the method computes variance components for gene-by-gene
#' interaction and genetic relatedness using the MQS estimator. The process is
#' parallelized using OpenMP when `n_threads > 1`.
#'
#' The memory requirements and computation time scaling can be optimized through
#' the parameters `chunk_size`, `n_randvecs`, and `n_blocks`.
#'
#' **Mask Format Requirements**
#'
#' The mask file format is an HDF5 file used for storing index data for
#' the masking process. This format supports data retrieval by index.
#' Below are the required groups and datasets within the HDF5 file:
#'
#' The required group names can be configured as input parameters.
#' The defaults are described below.
#'
#' - **Groups**:
#'   - `ld`: Stores SNPs in LD with the focal SNP. These SNPs will be **excluded**.
#'   - `gxg`: Stores indices of SNPs that the marginal epistasis test is conditioned on. These SNPs will be **included**.
#'
#' - **Datasets**:
#'   - `ld/<j>`: For each focal SNP `<j>`, this dataset contains indices of SNPs
#'     in the same LD block as that SNP. These SNPs will be **excluded** from the gene-by-gene interaction covariance matrix.
#'   - `gxg/<j>`: For each focal SNP `<j>`, this dataset contains indices of SNPs to **include** in the
#'     the gene-by-gene interaction covariance matrix for focal SNP `<j>`.
#'
#' **Important**: All indices in the mask file data are **zero-based**, matching the zero-based indices of the PLINK `.bim` file.
#'
#' @examples
#' plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed", package="sme"))
#' pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="sme")
#' mask_file <- ""
#'
#' # Parameter inputs
#' chunk_size <- 10
#' n_randvecs <- 10
#' n_blocks <- 10
#' n_threads <- 1
#'
#' # 1-based Indices of SNPs to be analyzed
#' n_snps <- 100
#' snp_indices <- 1:n_snps
#'
#' sme_result <- sme(
#'   plink_file,
#'   pheno_file,
#'   mask_file,
#'   snp_indices,
#'   chunk_size,
#'   n_randvecs,
#'   n_blocks,
#'   n_threads
#' )
#' head(sme_result$summary)
#'
#' @useDynLib SME
#' @import Rcpp
#' @import RcppEigen
#' @import dplyr
#' @importFrom stats pnorm
#' @importFrom tidyr pivot_longer
#' @importFrom utils read.delim
#' @export
sme <-
  function(plink_file,
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
           log_level = "WARNING") {
    logging::logReset()
    logging::basicConfig(level = log_level)
    log <- logging::getLogger("sme")

    n_gxg_indices <- length(gxg_indices)
    log$debug("Number of gxg indices: %d", n_gxg_indices)

    bim_file <- paste0(plink_file, ".bim")
    fam_file <- paste0(plink_file, ".fam")
    n_snps <- count_snps_bim(bim_file)
    n_samples <- count_samples(pheno_file)
    n_fam_lines <- count_fam(fam_file)

    log$debug("Dataset: %s", plink_file)
    log$debug("Trait file: %s", pheno_file)
    log$debug("Mask file: %s", mask_file)
    log$debug("Number of samples: %d", n_samples)
    log$debug("Number of SNPs: %d", n_snps)
    log$debug("Number of random vectors: %d", n_randvecs)
    log$debug("Number of blocks: %d", n_blocks)

    if (check_openmp()) {
      log$info("openMP is enabled")
      log$info("Number of requested threads: %d", n_threads)
    }

    if (n_samples != n_fam_lines) {
      stop("Number of samples in fam file and pheno file do not match.")
    }

    mem_req <- approximate_memory_requirements(n_samples, n_snps, n_blocks, n_randvecs, chunk_size)
    log$debug("Estimated memory requirement: %.2f GB", mem_req)

    if (is.null(gxg_indices)) {
      gxg_indices <- c(1:n_snps)
    }

    if (is.null(chunk_size)) {
      n_chunks <- ceiling(n_gxg_indices / n_threads)
      log$debug("No chunk size specified. Using %d chunks.", n_chunks)
    } else {
      n_chunks <- ceiling(n_gxg_indices / chunk_size)
      log$debug("Chunk size set to %d. Using %d chunks.", chunk_size, n_chunks)
    }

    if (n_chunks > 1) {
      chunks <- split(gxg_indices,
                      cut(seq_along(gxg_indices), n_chunks, labels = FALSE))
    } else {
      chunks <- list(`1` = gxg_indices)
    }

    VC <- NULL
    SE <- NULL
    TIME <- NULL

    for (i in seq_along(chunks)) {
      chunk <- chunks[[i]]
      result <-
        sme_cpp(
          plink_file,
          pheno_file,
          mask_file,
          n_randvecs,
          n_blocks,
          rand_seed,
          chunk - 1,
          # R is 1-indexed, C++ is 0-indexed
          n_threads,
          gxg_h5_group,
          ld_h5_group
        )
      VC <- rbind(VC, result$vc_estimate)
      SE <- rbind(SE, result$vc_se)
      TIME <- c(TIME, result$duration)
    }
    total_duration <- sum(TIME)
    average_duration <- total_duration / n_gxg_indices
    log$debug("Total computation time: %f seconds",
                  total_duration)
    log$debug("Average computation time per SNP: %f seconds",
              average_duration)

    z_score <- (VC / SE)
    p_values <- (1 - pnorm(z_score)) # one sided test

    pve <- VC / apply(VC, 1, sum)

    bim <- read.delim(bim_file, header = FALSE)
    colnames(bim) <- c("chromosome", "variant_id", "genetic_distance", "position", "allele1", "allele2")

    id <- bim$variant_id[gxg_indices]
    chromosome <- bim$chromosome[gxg_indices]
    position <- bim$position[gxg_indices]
    vc_names <- c("id", "grm", "gxg", "error")
    component_col <- "component"
    summary <-
      data.frame(
        id = id,
        index = gxg_indices,
        chromosome = chromosome,
        position = position,
        p = p_values[, 2],
        pve = pve[, 2],
        vc = VC[, 2],
        se = SE[, 2]
      )
    pve <- cbind(id, as.data.frame(pve))
    p_values <- cbind(id, as.data.frame(p_values))
    vc <- cbind(id, as.data.frame(VC))
    se <- cbind(id, as.data.frame(SE))
    colnames(pve) <- vc_names
    colnames(p_values) <- vc_names
    colnames(vc) <- vc_names
    colnames(se) <- vc_names
    result$p <- as_tibble(p_values)
    result$pve <- pivot_output(pve, component_col, "pve", vc_names[2:4])
    result$vc_estimate <- pivot_output(vc, component_col, "vc_estimate", vc_names[2:4])
    result$vc_se <- pivot_output(se, component_col, "vc_se", vc_names[2:4])
    result$summary <- as_tibble(summary)
    result$average_duration <- average_duration
    return(result)
  }


pivot_output <- function(df, names_to, values_to, cols) {
  as_tibble(df) %>% pivot_longer(cols = all_of(cols),
                                 names_to = names_to,
                                 values_to = values_to)
}
