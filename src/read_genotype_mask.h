#pragma once
#include <highfive/H5Easy.hpp>
#include <string>

#include "mme.h"

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const std::string &gxg_h5_dataset,
                        const std::string &ld_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps);

void set_mask_values(const std::string &genotype_mask_file, int n_snps,
                     int gxg_i, const std::string &h5_group,
                     MatrixXdr &genotype_mask, int &mask_value);

bool check_group_exists(const std::string &file_path, const std::string &group);
