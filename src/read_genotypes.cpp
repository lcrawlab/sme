#include "read_genotypes.h"

void read_focal_snp(const string &filename, MatrixXdr &focal_genotype,
                    const int &focal_snp_index, const int &n_samples,
                    const int &n_snps, int &global_snp_index) {
  ifstream ifs(filename.c_str(), ios::in | ios::binary);
  char magic[3];
  metaData metadata = set_metadata(n_samples);
  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];

  if (global_snp_index < 0) {
    binary_read(ifs, magic);
  }

  int y[4];

  for (int i = 0; i <= focal_snp_index; i++) {
    global_snp_index++;
    ifs.read(reinterpret_cast<char *>(gtype),
             metadata.ncol * sizeof(unsigned char));
    // skip to the n_encoded of interest
    if (i == (focal_snp_index)) {
      for (int k = 0; k < metadata.ncol; k++) {
        unsigned char c = gtype[k];
        unsigned char mask = metadata.mask;
        extract_plink_genotypes(y, c, mask);
        int j0 = k * metadata.unitsperword;
        int ncol = metadata.ncol;
        int lmax = get_sample_block_size(n_samples, k, ncol);
        for (int l = 0; l < lmax; l++) {
          int j = j0 + l;
          int val = encoding_to_allelecount(y[l]);
          // set missing genotype to zero to prevent the software from freezing
          val = (val == -1) ? 0 : val;
          focal_genotype(j, 0) = val;
        }
      }
    }
  }

  delete[] gtype;
}

void normalize_genotype(MatrixXdr &focal_genotype, const int &n_samples) {
  double mean_sel_snp = focal_genotype.array().sum() / n_samples;
  double sd_sel_snp = sqrt((mean_sel_snp * (1 - (0.5 * mean_sel_snp))));
  focal_genotype.array() = focal_genotype.array() - mean_sel_snp;
  focal_genotype.array() = focal_genotype.array() / sd_sel_snp;
}

int get_sample_block_size(const int &n_samples, const int &k, const int &ncol) {
  // Handle number of individuals not being a multiple of 4
  int lmax = 4;
  if (k == ncol - 1) {
    lmax = n_samples % 4;
    lmax = (lmax == 0) ? 4 : lmax;
  }
  return lmax;
}

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value) {
  return stream.read(reinterpret_cast<char *>(&value), sizeof(T));
}

void extract_plink_genotypes(int *y, const unsigned char &c,
                             const unsigned char &mask) {
  // Extract PLINK genotypes
  y[0] = (c)&mask;
  y[1] = (c >> 2) & mask;
  y[2] = (c >> 4) & mask;
  y[3] = (c >> 6) & mask;
}

void read_genotype_block(std::istream &ifs, const int &block_size,
                         genotype &genotype_block, const int &n_samples,
                         int &global_snp_index, const metaData &metadata) {

  for (int i = 0; i < block_size; i++) {
    MatrixXdr genotype_matrix = MatrixXdr::Zero(n_samples, 1);
    read_snp(ifs, global_snp_index, genotype_matrix);
    encode_snp(genotype_block, genotype_matrix);
  }
}

void encode_snp(genotype &genotype_block, const MatrixXdr &genotype_matrix) {
  int n_samples = genotype_matrix.rows();
  for (int j = 0; j < n_samples; j++) {
    encode_genotypes(genotype_block, j, genotype_matrix(j, 0));
  }
  genotype_block.n_encoded++;
}
void read_snp(std::istream &ifs, int &global_snp_index,
              MatrixXdr &genotype_matrix) {
  int n_samples = genotype_matrix.rows();
  metaData metadata = set_metadata(n_samples);
  char magic[3];
  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];
  if (global_snp_index < 0) {
    binary_read(ifs, magic);
  }
  int y[4];
  global_snp_index++;
  ifs.read(reinterpret_cast<char *>(gtype),
           metadata.ncol * sizeof(unsigned char));

  for (int k = 0; k < metadata.ncol; k++) {
    unsigned char c = gtype[k];
    unsigned char mask = metadata.mask;
    extract_plink_genotypes(y, c, mask);
    int j0 = k * metadata.unitsperword;
    int ncol = metadata.ncol;
    int lmax = get_sample_block_size(n_samples, k, ncol);
    for (int l = 0; l < lmax; l++) {
      int j = j0 + l;
      int val = encoding_to_allelecount(y[l]);
      // set missing genotype to major allele
      val = (val == -1) ? 0 : val;
      genotype_matrix(j, 0) = val;
    }
  }
  delete[] gtype;
}

void skip_snp(std::istream &ifs, int &global_snp_index, int &n_samples) {
  metaData metadata = set_metadata(n_samples);
  char magic[3];
  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];
  if (global_snp_index < 0) {
    binary_read(ifs, magic);
  }
  global_snp_index++;
  ifs.read(reinterpret_cast<char *>(gtype),
           metadata.ncol * sizeof(unsigned char));
  delete[] gtype;
}

void encode_genotypes(genotype &genotype_block, int j, int val) {
  int snp_index;
  snp_index = genotype_block.n_encoded;
  int horiz_seg_no = snp_index / genotype_block.segment_size_hori;
  genotype_block.p[horiz_seg_no][j] =
      3 * genotype_block.p[horiz_seg_no][j] + val;
  // computing sum for every snp to compute mean
  genotype_block.columnsum[snp_index] += val;
}

int encoding_to_allelecount(const int &value) {
  // Extract  PLINK coded genotype and convert into 0/1/2
  // PLINK coding:
  // 00->0
  // 01->missing
  // 10->1
  // 11->2
  std::string missing_message =
      "Missing genotype in PLINK file. This method requires fully genotyped "
      "data. Please impute missing genotypes.";
  switch (value) {
  case 0:
    return 0;
  case 1:
    // Rcpp::warning(missing_message); // Rcpp::warning is not thread safe.
    // TODO: remove?
    return -1;
  case 2:
    return 1;
  case 3:
    return 2;
  default:
    // Handle invalid input
    // Rcpp::warning(missing_message);
    return -1; // or any other default value you prefer
  }
}
