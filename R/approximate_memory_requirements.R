#' Estimate Memory Requirements for MME Routine
#'
#' This function provides an approximate estimate of the memory requirements 
#' (in gigabytes) for running the Multimodal Marginal Epistasis (MME) routine 
#' based on input parameters such as the number of samples, SNPs, and other configurations.
#'
#' @param n_samples Integer. The number of samples in the dataset.
#' @param n_snps Integer. The total number of SNPs in the dataset.
#' @param n_blocks Integer. The number of genotype blocks used to partition SNPs. 
#'   Affects the size of encoded genotype segments.
#' @param n_randvecs Integer. The number of random vectors used for stochastic 
#'   trace estimation. Affects memory for operations involving random vectors.
#' @param chunksize Integer. The number of focal SNPs processed per chunk.
#'
#' @return Numeric. The approximate memory requirement (in gigabytes) for the
#' MME routine.
#'
#' @details
#' The function calculates memory usage by summing the contributions from 
#' various components used in the MME routine, including:
#' - Variance component estimates (`vc_estimates`)
#' - Phenotype-related matrices
#' - Random vector-based computations
#' - Genotype objects and block statistics
#' - Gene-by-gene interaction masks
#'
#' The estimated memory requirement is derived from the data dimensions 
#' and operational needs, and it provides a guideline for configuring resources 
#' for the analysis.
#'
#' @examples
#' n_samples <- 1e5
#' n_snps <- 1e6
#' n_blocks <- 100
#' n_randvecs <- 100
#' chunksize <- 10
#' approximate_memory_requirements(n_samples, 
#'                                 n_snps, 
#'                                 n_blocks, 
#'                                 n_randvecs, 
#'                                 chunksize)
#'
#' @export
approximate_memory_requirements <- function(n_samples,
                                            n_snps,
                                            n_blocks,
                                            n_randvecs,
                                            chunksize) {
  n_encoded <- ceiling(n_snps / n_blocks)

  # VC - Matrix: (n_gxg_idx, n_variance_components + 1)
  # SE - Matrix: (n_gxg_idx, n_variance_components + 1)
  vc_estimates <- 2 * n_snps * 3 # point estimate and se for each component

  # pheno_mask - Matrix: (n_samples, 1)
  # pheno - Matrix: (n_samples, 1)
  # gxg_pheno - Matrix: (n_samples, 1)
  # snp_matrix - Matrix: (n_samples, 1)
  # focal_snp_gtype - Matrix: (n_samples, 1)
  # collect_XXy - Matrix: (n_samples, 1)
  # collect_Gy - Matrix: (n_samples, n_gxg_idx)
  # focal_snps_matrix - Matrix: (n_samples, n_gxg_idx)
  # collect_XXUy - Matrix: (n_samples, 
  #       (n_variance_components + 1) * (n_variance_components + 1) * n_gxg_idx)
  phenotype_like <- n_samples * (6 + 2 * chunksize + 9 * chunksize)

  # random_vectors - Matrix: (n_samples, n_randvecs)
  # gxg_random_vectors - Matrix: (n_samples, n_randvecs)
  # temp_grm - Matrix: (n_samples, n_randvecs)
  # temp_gxg - Matrix: (n_samples, n_randvecs)
  # XXz - Matrix: (n_samples, n_randvecs)
  # GxGz - Matrix: (n_samples, n_randvecs * n_gxg_idx)
  randomvec_like <- n_samples * n_randvecs * (5 + chunksize)

  # grm_genotype_block - genotype object
  # gxg_genotype_blocks - Vector of genotype objects: (n_gxg_idx)
  segment_size_hori <- floor(log(n_samples) / log(3)) - 2
  n_segments_hori <- ceiling(n_encoded / segment_size_hori)
  block_stats <- 2 * n_encoded * (chunksize + 1) # *2 for mean and variance

  # 4 bytes for int;
  gt_objects <- (n_segments_hori * n_samples) * (chunksize + 1) / 2

  # binary_gxg_mask - Matrix: (n_snps, n_gxg_idx)
  mask <- n_snps * chunksize

  # yXXy - Matrix: (1, 1)
  # yGxGy - Matrix: (n_gxg_idx, 1)
  # block_sizes - Vector<int>: (n_blocks)
  # n_gxg_snps_list - Vector<int>: (n_gxg_idx)
  # point_est - Matrix: (n_variance_components + 1, 1)
  # q - Matrix: (n_variance_components + 1, 1)
  # S - Matrix: (n_variance_components + 1, n_variance_components + 1)
  # cov_q - Matrix: (n_variance_components + 1, n_variance_components + 1)
  # invS - Matrix: (n_variance_components + 1, n_variance_components + 1)
  # cov_sigma - Matrix: (n_variance_components + 1, n_variance_components + 1)
  #
  # these are so small they can be neglected

  total <-  vc_estimates + 
            phenotype_like + 
            randomvec_like + 
            gt_objects + 
            block_stats + 
            mask
  return(total * 8 / 1024 / 1024 / 1024)
}
