#include "../../../mmer.Rcheck/00_pkg_src/mmer/src/mme.h"
#include <string>

int main() {
    std::string plink_file =
            "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed";
    std::string pheno_file =
            "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed.phen";
  std::string mask_file = ""; // "/Users/jds/Downloads/test100k/hdf5_mask.h5"
  std::string gxg_h5_group = "gxg";
  std::string ld_h5_group = "ld";
  int n_randvecs = 10;
  int focal_snp_index = 1;
  int n_blocks = 100;
  int rand_seed = 123;
  int n_threads = 1;
    std::vector<int> index_vector = {
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10};
  Rcpp::List results = mme_cpp(plink_file, pheno_file,
                               mask_file, plink_file,
                                n_randvecs, n_blocks, rand_seed, index_vector, n_threads, gxg_h5_group, ld_h5_group);
  return 0;
}
// 5 threads:
// - 5 snps 93 sec 18 / snp
// - 10 snps 154 sec - 15 / snp
// - 15 snps 227 sec - 15 / snp