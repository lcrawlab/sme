library(mmer)
library(tidyr)

plink_file <- "/Users/jds/data/ukbb/c12_100k-samples_020k-snps"
pheno_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps.pheno"
# plink_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed"
# pheno_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed.phen"
covariate_file <- ""
mask_file <- ""
mask_files <- c("", "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_masksize500.h5")
log_level <- "DEBUG"
n_blocks <- 100
rand_seed <- 123

snp_indices <- 1:90

chunksize <- c(90)
n_threads <- c(10)
n_randvecs <- c(5)

parameter_grid <- crossing(chunksize = chunksize,
                              n_threads = n_threads,
                              n_randvecs = n_randvecs,
                              mask_file = mask_files)

duration <- parameter_grid
duration$average_duration <- rep(0, nrow(parameter_grid))


for (i in 1:nrow(parameter_grid)) {
  # i <- 2
  result <- mme(plink_file,
                 pheno_file,
                 covariate_file,
                 parameter_grid$mask_file[i],
                 snp_indices,
                 parameter_grid$chunksize[i],
                 parameter_grid$n_randvecs[i],
                 n_blocks,
                 parameter_grid$n_threads[i],
                 rand_seed,
                 log_level)
  duration$average_duration[i] <- result$average_duration
}

print(duration)

