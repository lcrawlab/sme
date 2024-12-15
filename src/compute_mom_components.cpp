//
// Created by Julian Stamp on 4/2/24.
//

#include "compute_mom_components.h"

void compute_mom_components(int n_randvecs, int n_variance_components,
                            const MatrixXdr &pheno,
                            const MatrixXdr &random_vectors, MatrixXdr &XXz,
                            MatrixXdr &GxGz, const double &yXXy,
                            const double &yGxGy,
                            const std::vector<int> &n_snps_variance_component,
                            int n_samples_mask, MatrixXdr &S, MatrixXdr &q) {
  std::vector<MatrixXdr *> random_trace_components;
  random_trace_components.push_back(&XXz);
  random_trace_components.push_back(&GxGz);
  MatrixXdr A_trs(n_variance_components, n_variance_components);
  MatrixXdr b_trk(n_variance_components, 1);
  MatrixXdr c_yky(n_variance_components, 1);

  MatrixXdr yVy = MatrixXdr::Zero(2, 1);
  yVy(0, 0) = yXXy;
  yVy(1, 0) = yGxGy;

  for (int i = 0; i < n_variance_components; i++) {
    c_yky(i, 0) = yVy(i, 0) / n_snps_variance_component[i];

    for (int j = i; j < n_variance_components; j++) {
      MatrixXdr B1 = *random_trace_components[i];
      MatrixXdr B2 = *random_trace_components[j];

      double trkij = (B1.array() * B2.array()).colwise().sum().sum() /
                     n_snps_variance_component[i] /
                     n_snps_variance_component[j] / n_randvecs;
      A_trs(i, j) = A_trs(j, i) = trkij;
    }
  }

  MatrixXdr B1 = random_vectors.array() * GxGz.array();
  b_trk << n_samples_mask,
      B1.sum() / (n_snps_variance_component[1] * n_randvecs);
  S << A_trs, b_trk, b_trk.transpose(), n_samples_mask;
  q << c_yky, (pheno.array() * pheno.array()).sum();
}