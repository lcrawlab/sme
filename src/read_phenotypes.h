#include "mme.h"
#include <RcppEigen.h>
#include <fstream>

using namespace std;

void read_phenotypes(const int &n_samples, const std::string &filename,
                     MatrixXdr &pheno, MatrixXdr &mask);
