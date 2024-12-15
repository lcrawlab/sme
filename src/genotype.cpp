/*
 * Substantial parts of this file were published under MIT license by other
 * authors. Find the original license notice below.
 */

/* MIT License
 *
 * Copyright (c) 2024 sriramlab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "genotype.h"

using namespace std;

double genotype::get_col_mean(int snpindex) const {
  double temp = allelecount_means(snpindex, 0); // replaced above return
  // but did not seem to change the results
  return temp;
}

double genotype::get_col_std(int snpindex) const {
  double p_i = get_col_mean(snpindex);
  if (p_i == 0 || p_i == 2)
    return 1.0;
  double temp = sqrt(p_i * (1 - (0.5 * p_i)));
  return temp;
}

void genotype::clear_block() {
  std::vector<std::vector<int>>().swap(p);
  n_encoded = 0;
  columnsum.resize(1, 0);
  allelecount_means.resize(1, 1);
  allelecount_stds.resize(1, 1);
}

void genotype::set_block_parameters(const int &n, const int &b) {
  segment_size_hori = floor(log(n) / log(3)) - 2; // object of the mailman
  // TODO: the above line introduces a minimum size limit. Investigate the
  //  limitations
  n_segments_hori = ceil(b * 1.0 / (segment_size_hori * 1.0));
  p.resize(n_segments_hori, std::vector<int>(n, 0));
  n_encoded = 0; // will be updated in encode
  n_snps = b;
  n_samples = n;

  columnsum.resize(b, 1);
  for (int index_temp = 0; index_temp < b; index_temp++) {
    columnsum[index_temp] = 0;
  }
}

void genotype::compute_block_stats() {
  allelecount_stds.resize(n_encoded, 1);
  allelecount_means.resize(n_encoded, 1);
  for (int i = 0; i < n_encoded; i++) {
    allelecount_means(i, 0) = (double)columnsum[i] / n_samples;
    allelecount_stds(i, 0) =
        1 /
        sqrt((allelecount_means(i, 0) * (1 - (0.5 * allelecount_means(i, 0)))));
  }
}

void genotype::encode_snp(const MatrixXdr &snp_matrix) {
  int horiz_seg_no = n_encoded / segment_size_hori;
  for (int j = 0; j < n_samples; j++) {
    int val = snp_matrix(j, 0);
    p[horiz_seg_no][j] = 3 * p[horiz_seg_no][j] + val;
    // computing sum for every snp to compute mean
    columnsum[n_encoded] += val;
  }
  n_encoded++;
}
