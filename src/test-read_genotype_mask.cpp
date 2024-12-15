/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>

#include "testing_utils.h"

context("C++ test reading genotype mask file") {
  test_that("Mask sets correct snps to 1 and 0") {
    correctTestFiles(test_csv, test_bed, test_pheno, test_h5, test_ld_h5);
    // given
    std::string genotype_mask_file;
    int n_snps = 100;
    int gxg_i = 5;
    std::string gxg_h5_dataset = "gxg";
    std::string ld_h5_dataset = "ld";
    MatrixXdr genotype_mask = MatrixXdr::Zero(n_snps, 1);
    int n_gxg_snps;
    // when
    read_genotype_mask(test_h5, n_snps, gxg_i, gxg_h5_dataset, ld_h5_dataset,
                       genotype_mask, n_gxg_snps);

    // then
    expect_true(genotype_mask.sum() == n_snps - 1);
    expect_true(genotype_mask(gxg_i, 0) == 0);
  }

  test_that("Mask sets correct snps to 1 and 0") {
    // given
    std::string test_h5 = "";
    int n_snps = 100;
    int gxg_i = 5;
    std::string gxg_h5_dataset = "gxg";
    std::string ld_h5_dataset = "ld";
    MatrixXdr genotype_mask = MatrixXdr::Zero(n_snps, 1);
    int n_gxg_snps;
    // when
    read_genotype_mask(test_h5, n_snps, gxg_i, gxg_h5_dataset, ld_h5_dataset,
                       genotype_mask, n_gxg_snps);

    // then
    expect_true(genotype_mask.sum() == n_snps);
    expect_true(genotype_mask(gxg_i, 0) == 1);
  }

  test_that("LD mask is incorporated") {
    correctTestFiles(test_csv, test_bed, test_pheno, test_h5, test_ld_h5);
    // given
    std::string genotype_mask_file;
    int n_snps = 100;
    int gxg_i = 5;
    std::string gxg_h5_dataset = "gxg";
    std::string ld_h5_dataset = "ld";
    MatrixXdr genotype_mask = MatrixXdr::Zero(n_snps, 1);
    int n_gxg_snps;
    // when
    read_genotype_mask(test_ld_h5, n_snps, gxg_i, gxg_h5_dataset, ld_h5_dataset,
                       genotype_mask, n_gxg_snps);
    // then
    expect_true(genotype_mask.sum() == n_snps - 1 - 10);
    expect_true(genotype_mask(gxg_i, 0) == 0);
  }
}
