#include "read_genotype_mask.h"

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const std::string &gxg_h5_dataset,
                        const std::string &ld_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps) {
  int include = 1;
  int exclude = 0;

  if (genotype_mask_file != "") {
    bool has_gxg_data = check_group_exists(genotype_mask_file, gxg_h5_dataset);
    bool has_ld_data = check_group_exists(genotype_mask_file, ld_h5_dataset);

    if (has_gxg_data) {
      set_mask_values(genotype_mask_file, n_snps, gxg_i, gxg_h5_dataset,
                      genotype_mask, include);
    } else {
      genotype_mask = MatrixXdr::Ones(n_snps, 1);
      n_gxg_snps = n_snps - 1;
    }

    if (has_ld_data) {
      set_mask_values(genotype_mask_file, n_snps, gxg_i, ld_h5_dataset,
                      genotype_mask, exclude);
    }

    n_gxg_snps = genotype_mask.sum();
  } else {
    genotype_mask = MatrixXdr::Ones(n_snps, 1);
    n_gxg_snps = n_snps - 1;
  }
}

void set_mask_values(const std::string &genotype_mask_file, int n_snps,
                     int gxg_i, const std::string &h5_group,
                     MatrixXdr &genotype_mask, int &mask_value) {
  std::vector<int> genotype_mask_indices;
  H5Easy::File file(genotype_mask_file, H5Easy::File::ReadOnly);
  genotype_mask_indices = H5Easy::load<std::vector<int>>(
      file, h5_group + "/" + std::to_string(gxg_i));
  for (int index : genotype_mask_indices) {
    if (index >= 0 && index < n_snps) {
      genotype_mask(index, 0) = mask_value;
    }
  }
}

/**
 * @brief Check for the existence of "ld" group in an HDF5 file
 *
 * @param file_path The path to the HDF5 file
 * @return true if the group exists in the data set, false otherwise
 */
bool check_group_exists(const std::string &file_path,
                        const std::string &group) {
  try {
    // Open the HDF5 file
    HighFive::File file(file_path, HighFive::File::ReadOnly);
    return file.exist(group);
  } catch (const HighFive::Exception &e) {
    return false;
  }
}
