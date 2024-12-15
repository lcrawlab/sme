#ifndef MMER_TESTING_UTILS_H
#define MMER_TESTING_UTILS_H

#endif // MMER_TESTING_UTILS_H

#include "allocate_memory.h"
#include "computation.h"
#include "compute_covariance_q.h"
#include "compute_mom_components.h"
#include "genotype.h"
#include "initialize_random_vectors.h"
#include "mme.h"
#include "read_genotype_mask.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include "set_metadata.h"

#include <fstream>
#include <iostream>
#include <unistd.h>

// this should come from one configuration file
extern std::string testdata_dir;
extern std::string checkdata_dir;
extern std::string test_bed;
extern std::string test_csv;
extern std::string test_pheno;
extern std::string test_h5;
extern std::string test_ld_h5;
extern std::string test_bed_2;

extern bool is_test;

extern double tolerance;
extern int n_samples;
extern int block_size;
extern metaData metadata;
extern metaData metadata2;

MatrixXdr readCSVToMatrixXdr(const std::string &filename);

bool fileExists(const std::string &path);

void correctTestFiles(std::string &test_csv, std::string &test_bed,
                      std::string &test_pheno, std::string &test_h5,
                      std::string &test_ld_h5);
