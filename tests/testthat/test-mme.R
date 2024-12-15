test_that("mme end-to-end no mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed",
                                               package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- ""
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "WARNING"

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947,
      0.424659762, -0.034297866, 0.586712315,
      0.423531782, -0.057827236, 0.610711311,
      0.424752022, -0.040029493, 0.589124766,
      0.422525234, -0.038636283, 0.590873964,
      0.423146054, 0.022630877, 0.525568865,
      0.424032481, -0.007207079, 0.557439240,
      0.423783848, 0.097367174, 0.448923311,
      0.424160268, -0.014771179, 0.565349071,
      0.423546323, 0.017600500, 0.535181598
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- mme(plink_file,
                 pheno_file,
                 mask_file,
                 snp_indices,
                 chunksize,
                 n_randvecs,
                 n_blocks,
                 n_threads,
                 gxg_h5_group,
                 ld_h5_group,
                 rand_seed,
                 log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end with mask", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed",
                                               package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- system.file("testdata", "test.h5", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 1

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.414218, 0.232245, 0.363624,
      0.425354, -0.0394186, 0.594756,
      0.42462, -0.0568018, 0.607989,
      0.42503, -0.0384796, 0.587712,
      0.424906, 0.0268949, 0.521808,
      0.415192, 0.0974395, 0.43878,
      0.424643, -0.0158059, 0.565519,
      0.42001, 0.146484, 0.390797,
      0.424789, -0.0199609, 0.570214,
      0.423682, 0.0074807, 0.5438
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )

  id <- sprintf("rs%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- mme(plink_file,
                 pheno_file,
                 mask_file,
                 snp_indices,
                 chunksize,
                 n_randvecs,
                 n_blocks,
                 n_threads,
                 gxg_h5_group,
                 ld_h5_group,
                 rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-1)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end no mask only one gxg idx", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed",
                                               package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  mask_file <- ""
  chunksize <- 3
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "WARNING"

  snp_indices <- c(3)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed,
                log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end no mask - chunksize 1", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed",
                                               package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- ""
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 1
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 3
  log_level <- "WARNING"

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.422807472, 0.030770829, 0.525349947,
      0.424659762, -0.034297866, 0.586712315,
      0.423531782, -0.057827236, 0.610711311,
      0.424752022, -0.040029493, 0.589124766,
      0.422525234, -0.038636283, 0.590873964,
      0.423146054, 0.022630877, 0.525568865,
      0.424032481, -0.007207079, 0.557439240,
      0.423783848, 0.097367174, 0.448923311,
      0.424160268, -0.014771179, 0.565349071,
      0.423546323, 0.017600500, 0.535181598
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )
  id <- sprintf("rs%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed,
                log_level)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-2)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})

test_that("mme end-to-end but with mask - chunksize 1", {
  # given
  plink_file <- gsub("\\.bed", "", system.file("testdata", "test.bed",
                                               package="mmer"))
  pheno_file <- system.file("testdata", "test_h2_0.5.pheno", package="mmer")
  mask_file <- system.file("testdata", "test.h5", package="mmer")
  gxg_h5_group <- "gxg"
  ld_h5_group <- "ld"
  chunksize <- 1
  n_randvecs <- 10
  n_blocks <- 10
  rand_seed <- 123
  n_threads <- 1

  snp_indices <- c(3, 8, 9, 13, 16, 19, 29, 34, 93, 97)

  expected_est <- matrix(
    c(
      0.414218, 0.232245, 0.363624,
      0.425354, -0.0394186, 0.594756,
      0.42462, -0.0568018, 0.607989,
      0.42503, -0.0384796, 0.587712,
      0.424906, 0.0268949, 0.521808,
      0.415192, 0.0974395, 0.43878,
      0.424643, -0.0158059, 0.565519,
      0.42001, 0.146484, 0.390797,
      0.424789, -0.0199609, 0.570214,
      0.423682, 0.0074807, 0.5438
    ), ncol = 3, byrow = TRUE
  )
  expected_se <- matrix(
    c(
      0.12076409, 0.07483502, 0.10296534,
      0.12092711, 0.02880959, 0.08871578,
      0.12084017, 0.04131487, 0.10085628,
      0.12135710, 0.04443874, 0.09379178,
      0.12059002, 0.04714957, 0.09737975,
      0.12084237, 0.04549436, 0.09319999,
      0.12096347, 0.05326928, 0.09342799,
      0.12158889, 0.10457918, 0.12324774,
      0.12089678, 0.04467541, 0.09617960,
      0.12099050, 0.07292024, 0.10320580
    ), ncol = 3, byrow = TRUE
  )

  id <- sprintf("rs%d", snp_indices)
  vc_names <- c("id", "grm", "gxg", "error")
  vc_df <- cbind(id, as.data.frame(expected_est))
  colnames(vc_df) <- vc_names
  vc_df <- pivot_output(vc_df,
                        "component",
                        "vc_estimate",
                        vc_names[2:4])
  se_df <- cbind(id, as.data.frame(expected_se))
  colnames(se_df) <- vc_names
  se_df <- pivot_output(se_df,
                        "component",
                        "vc_se",
                        vc_names[2:4])
  # when
  result <- mme(plink_file,
                pheno_file,
                mask_file,
                snp_indices,
                chunksize,
                n_randvecs,
                n_blocks,
                n_threads,
                gxg_h5_group,
                ld_h5_group,
                rand_seed)
  observed_est <- result$vc_estimate
  observed_se <- result$vc_se

  # then
  expect_equal(observed_est, vc_df, tolerance = 1e-1)
  expect_equal(observed_se, se_df, tolerance = 1e-1)
})
