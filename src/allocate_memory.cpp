//
// Created by Julian Stamp on 3/29/24.
//

#include "allocate_memory.h"

void allocate_memory(int n_randvecs, const genotype &genotype_block,
                     double *&partialsums, double *&sum_op, double *&yint_e,
                     double *&yint_m, double **&y_e, double **&y_m) {
  int hsegsize = genotype_block.segment_size_hori;
  int hsize = pow(3, hsegsize); // = log_3(n)

  partialsums =
      new double[n_randvecs](); // value initialization with () to value 0
  sum_op = new double[n_randvecs]();
  yint_e = new double[hsize * n_randvecs];
  yint_m = new double[hsize * n_randvecs];
  memset(yint_m, 0, hsize * n_randvecs * sizeof(double));
  memset(yint_e, 0, hsize * n_randvecs * sizeof(double));

  y_e = new double *[genotype_block.n_samples];
  for (int i = 0; i < genotype_block.n_samples; i++) {
    y_e[i] = new double[n_randvecs];
    memset(y_e[i], 0, n_randvecs * sizeof(double));
  }

  y_m = new double *[hsegsize];
  for (int i = 0; i < hsegsize; i++)
    y_m[i] = new double[n_randvecs]();
}

void deallocate_memory(double *partialsums, double *sum_op, double *yint_e,
                       double *yint_m, double **y_e, double **y_m,
                       const genotype &genotype_block) {
  int hsegsize = genotype_block.segment_size_hori;
  delete[] sum_op;
  delete[] partialsums;
  delete[] yint_e;
  delete[] yint_m;
  for (int i = 0; i < hsegsize; i++)
    delete[] y_m[i];
  delete[] y_m;

  for (int i = 0; i < genotype_block.n_samples; i++)
    delete[] y_e[i];
  delete[] y_e;
}