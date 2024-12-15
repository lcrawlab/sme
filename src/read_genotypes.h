#include <RcppEigen.h>
#include <fstream>
#include <iostream>

#include "mme.h"
#include "set_metadata.h"

using namespace std;

void read_focal_snp(const string &filename, MatrixXdr &focal_genotype,
                    const int &focal_snp_index, const int &n_samples,
                    const int &n_snps, int &global_snp_index);

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value);

void read_genotype_block(std::istream &ifs, const int &block_size,
                         genotype &genotype_block, const int &n_samples,
                         int &global_snp_index, const metaData &metadata);

void read_snp(std::istream &ifs, int &global_snp_index,
              MatrixXdr &genotype_matrix);

void skip_snp(std::istream &ifs, int &global_snp_index, int &n_samples);

void encode_snp(genotype &genotype_block, const MatrixXdr &genotype_matrix);

void extract_plink_genotypes(int *y, const unsigned char &c,
                             const unsigned char &mask);

int encoding_to_allelecount(const int &value);

void normalize_genotype(MatrixXdr &focal_genotype, const int &n_samples);

int get_sample_block_size(const int &n_samples, const int &k, const int &ncol);

void encode_genotypes(genotype &genotype_block, int j, int val);
