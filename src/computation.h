#pragma once
#include "mme.h"

MatrixXdr compute_XXz(const MatrixXdr &Z_b, const MatrixXdr &phenotype_mask,
                      const int &n_randvecs, const genotype &genotype_block);

double compute_yXXy(genotype &genotype_block, const MatrixXdr &y_vec);

void multiply_y_pre_fast(const genotype &g, const MatrixXdr &op, MatrixXdr &res,
                         const bool &subtract_means, double *&sum_op,
                         double *&yint_m, double **&y_m, double *&partialsums);

void multiply_y_post_fast(const genotype g, const MatrixXdr &op, MatrixXdr &res,
                          bool subtract_means, double *&yint_e, double **&y_e,
                          int n);
