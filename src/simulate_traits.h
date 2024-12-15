// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "boost/random.hpp"
#include "count_data.h"
#include "genotype.h"
#include "read_genotypes.h"
#ifndef MMER_SIMULATE_TRAITS_H
#define MMER_SIMULATE_TRAITS_H

#endif // MMER_SIMULATE_TRAITS_H

std::vector<int> draw_random_ints(std::vector<int> numbers, int x);

MatrixXdr draw_normal_effects(int n);

MatrixXdr &scale_component(float target_variance, MatrixXdr &component);

Rcpp::List simulate_traits(std::string plink_file, float additive_heritability,
                           float gxg_heritability,
                           std::vector<int> additive_snps,
                           std::vector<int> gxg_group_1,
                           std::vector<int> gxg_group_2);
