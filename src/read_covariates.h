#include "mme.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

int read_covariates(const bool &std, const int &Nind,
                    const std::string &filename, MatrixXdr &covariate);
