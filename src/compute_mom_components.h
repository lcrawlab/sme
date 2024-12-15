//
// Created by Julian Stamp on 4/2/24.
//

#ifndef MMER_COMPUTE_MOM_COMPONENTS_H
#define MMER_COMPUTE_MOM_COMPONENTS_H
#include "mme.h"
#include <vector>

void compute_mom_components(int n_randvecs, int n_variance_components,
                            const MatrixXdr &pheno,
                            const MatrixXdr &random_vectors, MatrixXdr &XXz,
                            MatrixXdr &GxGz, const double &yXXy,
                            const double &yGxGy,
                            const std::vector<int> &n_snps_variance_component,
                            int n_samples_mask, MatrixXdr &S, MatrixXdr &q);

#endif // MMER_COMPUTE_MOM_COMPONENTS_H
