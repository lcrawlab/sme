#pragma once
#ifndef GENOTYPE_H
#define GENOTYPE_H
#include "mme.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <istream>

class genotype {

public:
  std::vector<int> columnsum;
  MatrixXdr allelecount_means;
  MatrixXdr allelecount_stds;

  int n_encoded;
  int n_snps;
  int n_samples;
  int n_segments_hori;
  int segment_size_hori;
  std::vector<std::vector<int>> p;

  // Constructor to initialize all members
  genotype()
      : columnsum(), allelecount_means(MatrixXdr::Zero(1, 1)),
        allelecount_stds(MatrixXdr::Zero(1, 1)), n_encoded(0), n_snps(0),
        n_samples(0), n_segments_hori(0), segment_size_hori(0), p() {}

  double get_col_mean(int snpindex) const;
  double get_col_std(int snpindex) const;
  void clear_block();
  void set_block_parameters(const int &n, const int &b);
  void compute_block_stats();
  void encode_snp(const MatrixXdr &snp_matrix);
};

#endif
