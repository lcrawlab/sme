set.seed(123)
library(dplyr)
library(genio)
library(mmer)

plink_file <- tempfile()
n_samples <- 5000
n_snps <- 6000
maf <- 0.05 + 0.45 * runif(n_snps)
random_minor_allele_counts   <- (runif(n_samples * n_snps) < maf) + (runif(n_samples * n_snps) < maf)
allele_count_matrix <- matrix(random_minor_allele_counts,
                              nrow = n_samples,
                              ncol = n_snps,
                              byrow = TRUE,
)

# Create .fam and .bim data frames
fam <- data.frame(
  fam = sprintf("F%d", 1:n_samples),       # Family ID
  id = sprintf("I%d", 1:n_samples),       # Individual ID
  pat = 0,                        # Paternal ID
  mat = 0,                        # Maternal ID
  sex = sample(1:2, n_samples, replace = TRUE),
  pheno = -9
)

bim <- data.frame(
  chr = 1,               # Chromosome number
  id = sprintf("rs%d", 1:n_snps),   # SNP name
  posg = 0,                         # Genetic distance (cM)
  pos = 1:n_snps,          # Base-pair position
  ref = sample(c("A", "C", "G", "T"), n_snps, replace = T),          # Minor allele
  alt = sample(c("A", "C", "G", "T"), n_snps, replace = T)         # Major allele
)

# Write to .bed, .fam, and .bim files
write_plink(
  file = plink_file,
  X = t(allele_count_matrix),
  fam = fam,
  bim = bim
)

pheno_file <- tempfile()
additive_heritability <- 0.3
gxg_heritability <- 0.25
n_snps_gxg_group <- 5
n_snps_additive <- 100
additive_snps <- sort(sample(1:n_snps, n_snps_additive, replace = F))
gxg_group_1 <- sort(sample(additive_snps, n_snps_gxg_group, replace = F))
gxg_group_2 <- sort(sample(setdiff(additive_snps, gxg_group_1), n_snps_gxg_group, replace = F))

simulate_traits(
  plink_file,
  pheno_file,
  additive_heritability,
  gxg_heritability,
  additive_snps,
  gxg_group_1,
  gxg_group_2
)

mask_file <- ""
gxg_h5_group <- "gxg"
ld_h5_group <- "ld"
chunk_size <- 10
n_randvecs <- 10
n_blocks <- 10
rand_seed <- 123
n_threads <- 5
log_level <- "DEBUG"

getting_started <- mme(
  plink_file,
  pheno_file,
  mask_file,
  additive_snps, # includes the gxg snps
  chunk_size,
  n_randvecs,
  n_blocks,
  n_threads,
  gxg_h5_group,
  ld_h5_group,
  rand_seed,
  log_level
)

getting_started$summary <- getting_started$summary %>%
  mutate(true_gxg_snp = (index %in% c(gxg_group_1, gxg_group_2)))

usethis::use_data(getting_started, overwrite = TRUE)
