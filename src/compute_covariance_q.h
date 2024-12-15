//
// Created by Julian Stamp on 4/2/24.
//

#ifndef MMER_COMPUTE_COVARIANCE_Q_H
#define MMER_COMPUTE_COVARIANCE_Q_H

#include "mme.h"

void compute_covariance_q(int n_variance_components,
                          const MatrixXdr &collect_UVy,
                          const MatrixXdr &point_est, MatrixXdr &cov_q);

#endif // MMER_COMPUTE_COVARIANCE_Q_H
