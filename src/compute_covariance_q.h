//
// Created by Julian Stamp on 4/2/24.
//

#ifndef SME_COMPUTE_COVARIANCE_Q_H
#define SME_COMPUTE_COVARIANCE_Q_H

#include "sme.h"

void compute_covariance_q(int n_variance_components,
                          const MatrixXdr &collect_UVy,
                          const MatrixXdr &point_est, MatrixXdr &cov_q);

#endif // SME_COMPUTE_COVARIANCE_Q_H
