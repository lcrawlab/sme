#include "mme.h"
#include <string>

int main() {
    std::string plink_file =
            "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed";
    std::string pheno_file =
            "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed.phen";
    std::string mask_file = "";
    std::string gxg_h5_group = "gxg";
    std::string ld_h5_group = "ld";
    int n_randvecs = 100;
    int focal_snp_index = 1;
    int n_blocks = 100;
    int rand_seed = 123;
    int n_threads = 10;
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
            10 //,
            //                                   11,
            //                                   12,
            //                                   13,
            //                                   14,
            //                                   15,
    };
    Rcpp::List results =
            mme_cpp(plink_file, pheno_file, mask_file, plink_file, n_randvecs, n_blocks,
                     rand_seed, index_vector, n_threads, gxg_h5_group, ld_h5_group);
    return 0;
}
// 5 threads:
// - 5 snps 241 sec - 48 / snp
// - 15 snps 700 sec - 35 / snp