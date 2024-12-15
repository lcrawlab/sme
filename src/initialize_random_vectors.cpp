//
// Created by Julian Stamp on 3/29/24.
//

#include "initialize_random_vectors.h"

MatrixXdr &initialize_random_vectors(int n_randvecs, int rand_seed,
                                     const MatrixXdr &pheno_mask,
                                     MatrixXdr &random_vectors, int n_samples) {
  // TODO: signature should either return matrix or pass it by reference...
  //  clean this up!
  random_vectors = MatrixXdr::Zero(n_samples, n_randvecs);

  boost::mt19937 seedr;
  seedr.seed(rand_seed < 0 ? std::time(0) : rand_seed);
  boost::normal_distribution<> dist(0, 1);
  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
      z_vec(seedr, dist);

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      random_vectors(j, i) = z_vec();

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      random_vectors(j, i) = random_vectors(j, i) * pheno_mask(j, 0);
  return random_vectors;
}