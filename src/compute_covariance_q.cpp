//
// Created by Julian Stamp on 4/2/24.
//

#include "compute_covariance_q.h"

void compute_covariance_q(int n_variance_components,
                          const MatrixXdr &collect_UVy,
                          const MatrixXdr &point_est, MatrixXdr &cov_q) {
  cov_q.resize(n_variance_components + 1, n_variance_components + 1);
  for (int k = 0; k < (n_variance_components + 1); k++) {
    for (int l = 0; l < (n_variance_components + 1); l++) {
      double total = 0;
      for (int t = 0; t < (n_variance_components + 1); t++) {

        MatrixXdr UVy = collect_UVy.col((t * (n_variance_components + 1)) + l);
        MatrixXdr Vy = collect_UVy.col(
            n_variance_components * (n_variance_components + 1) + k);
        double temp = (UVy.array() * Vy.array()).sum();
        total += point_est(t, 0) * temp;
      }
      cov_q(k, l) = 2 * total;
    }
  }
}