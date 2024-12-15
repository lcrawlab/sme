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

#include "computation.h"
#include "allocate_memory.h"
#include "mailman.h"

MatrixXdr compute_XXz(const MatrixXdr &Z_b, const MatrixXdr &phenotype_mask,
                      const int &n_randvecs, const genotype &genotype_block) {

  double *partialsums;
  double *sum_op;
  double *yint_e;
  double *yint_m;
  double **y_e;
  double **y_m;
  const MatrixXdr &means = genotype_block.allelecount_means;
  const MatrixXdr &stds = genotype_block.allelecount_stds;
  const int &n_samples = genotype_block.n_samples;

  allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                  yint_m, y_e, y_m);

  MatrixXdr res = MatrixXdr::Zero(genotype_block.n_encoded, n_randvecs);
  multiply_y_pre_fast(genotype_block, Z_b, res, false, sum_op, yint_m, y_m,
                      partialsums);
  MatrixXdr zb_sum = Z_b.colwise().sum();

  for (int j = 0; j < genotype_block.n_encoded; j++)
    for (int k = 0; k < n_randvecs; k++) {
      res(j, k) = res(j, k) * stds(j, 0);
    }

  MatrixXdr resid = MatrixXdr::Zero(genotype_block.n_encoded, n_randvecs);
  MatrixXdr inter = means.cwiseProduct(stds);
  resid = inter * zb_sum;
  MatrixXdr inter_zb = MatrixXdr::Zero(genotype_block.n_encoded, n_randvecs);
  inter_zb = res - resid;

  for (int j = 0; j < genotype_block.n_encoded; j++)
    for (int k = 0; k < n_randvecs; k++) {
      inter_zb(j, k) = inter_zb(j, k) * stds(j, 0);
    }

  MatrixXdr new_zb = inter_zb.transpose();
  MatrixXdr new_res = MatrixXdr::Zero(n_randvecs, n_samples);

  multiply_y_post_fast(genotype_block, new_zb, new_res, false, yint_e, y_e,
                       n_samples);
  MatrixXdr zb_scale_sum = new_zb * means;
  MatrixXdr new_resid = zb_scale_sum * MatrixXdr::Ones(1, n_samples);

  // new zb
  MatrixXdr temp = new_res - new_resid; // does it apply it to every column?

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      temp(i, j) = temp(i, j) * phenotype_mask(j, 0);

  deallocate_memory(partialsums, sum_op, yint_e, yint_m, y_e, y_m,
                    genotype_block);

  return temp.transpose();
}

double compute_yXXy(genotype &genotype_block, const MatrixXdr &y_vec) {
  double *partialsums;
  double *sum_op;
  double *yint_e;
  double *yint_m;
  double **y_e;
  double **y_m;
  int Nz = 1;
  const MatrixXdr &means = genotype_block.allelecount_means;
  const MatrixXdr &stds = genotype_block.allelecount_stds;

  allocate_memory(Nz, genotype_block, partialsums, sum_op, yint_e, yint_m, y_e,
                  y_m);

  MatrixXdr res = MatrixXdr::Zero(genotype_block.n_encoded, 1);

  multiply_y_pre_fast(genotype_block, y_vec, res, false, sum_op, yint_m, y_m,
                      partialsums);

  res = res.cwiseProduct(stds);
  MatrixXdr resid = MatrixXdr::Zero(genotype_block.n_encoded, 1);
  resid = means.cwiseProduct(stds);
  resid = resid * y_vec.sum();

  MatrixXdr Xy = MatrixXdr::Zero(genotype_block.n_encoded, 1);
  Xy = res - resid;

  double yXXy = (Xy.array() * Xy.array()).sum();

  deallocate_memory(partialsums, sum_op, yint_e, yint_m, y_e, y_m,
                    genotype_block);

  return yXXy;
}

void multiply_y_pre_fast(const genotype &g, const MatrixXdr &op, MatrixXdr &res,
                         const bool &subtract_means, double *&sum_op,
                         double *&yint_m, double **&y_m, double *&partialsums) {
  bool var_normalize = false;
  int Ncol_op = op.cols();

  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    sum_op[k_iter] = op.col(k_iter).sum();
  }

  for (int seg_iter = 0; seg_iter < g.n_segments_hori - 1; seg_iter++) {
    mailman::fastmultiply(g.segment_size_hori, g.n_samples, Ncol_op,
                          g.p[seg_iter], op, yint_m, partialsums, y_m);
    int p_base = seg_iter * g.segment_size_hori;
    for (int p_iter = p_base;
         (p_iter < p_base + g.segment_size_hori) && (p_iter < g.n_snps);
         p_iter++) {
      for (int k_iter = 0; k_iter < Ncol_op; k_iter++)
        res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
    }
  }

  int last_seg_size = (g.n_snps % g.segment_size_hori != 0)
                          ? g.n_snps % g.segment_size_hori
                          : g.segment_size_hori;
  mailman::fastmultiply(last_seg_size, g.n_samples, Ncol_op,
                        g.p[g.n_segments_hori - 1], op, yint_m, partialsums,
                        y_m);
  int p_base = (g.n_segments_hori - 1) * g.segment_size_hori;
  for (int p_iter = p_base;
       (p_iter < p_base + g.segment_size_hori) && (p_iter < g.n_snps);
       p_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++)
      res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
  }

  if (!subtract_means) {
    return;
  }

  for (int p_iter = 0; p_iter < g.n_encoded; p_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
      res(p_iter, k_iter) =
          res(p_iter, k_iter) - (g.get_col_mean(p_iter) * sum_op[k_iter]);
      if (var_normalize)
        res(p_iter, k_iter) = res(p_iter, k_iter) / (g.get_col_std(p_iter));
    }
  }
}

void multiply_y_post_fast(const genotype g, const MatrixXdr &op, MatrixXdr &res,
                          bool subtract_means, double *&yint_e, double **&y_e,
                          int n) {
  bool var_normalize = false;
  double *partialsums = new double[0];

  MatrixXdr op_t;
  op_t = op.transpose();

  int Ncol_op = op_t.cols();

  if (var_normalize && subtract_means) {
    for (int p_iter = 0; p_iter < g.n_encoded; p_iter++) {
      for (int k_iter = 0; k_iter < Ncol_op; k_iter++)
        op_t(p_iter, k_iter) = op_t(p_iter, k_iter) / (g.get_col_std(p_iter));
    }
  }

  int seg_iter;
  for (seg_iter = 0; seg_iter < g.n_segments_hori - 1; seg_iter++) {
    mailman::fastmultiply_pre(g.segment_size_hori, g.n_samples, Ncol_op,
                              seg_iter * g.segment_size_hori, g.p[seg_iter],
                              op_t, yint_e, partialsums, y_e);
  }
  int last_seg_size = (g.n_snps % g.segment_size_hori != 0)
                          ? g.n_snps % g.segment_size_hori
                          : g.segment_size_hori;
  mailman::fastmultiply_pre(last_seg_size, g.n_samples, Ncol_op,
                            seg_iter * g.segment_size_hori, g.p[seg_iter], op_t,
                            yint_e, partialsums, y_e);
  delete[] partialsums;
  for (int n_iter = 0; n_iter < n; n_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
      res(k_iter, n_iter) = y_e[n_iter][k_iter];
      y_e[n_iter][k_iter] = 0;
    }
  }

  if (!subtract_means)
    return;

  double *sums_elements = new double[Ncol_op];
  std::fill(sums_elements, sums_elements + Ncol_op, 0.0);

  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    double sum_to_calc = 0.0;
    for (int p_iter = 0; p_iter < g.n_encoded; p_iter++)
      sum_to_calc += g.get_col_mean(p_iter) * op_t(p_iter, k_iter);
    sums_elements[k_iter] = sum_to_calc;
  }
  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    for (int n_iter = 0; n_iter < n; n_iter++)
      res(k_iter, n_iter) = res(k_iter, n_iter) - sums_elements[k_iter];
  }
  delete[] sums_elements;
}
