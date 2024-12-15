#include <highfive/H5Easy.hpp>
#include <Rcpp.h>

// [[Rcpp::export]]
void createH5File(const std::string &filename) {
  HighFive::File file(filename, HighFive::File::Truncate);
}


// [[Rcpp::export]]
std::vector<int> readH5File(const std::string &filename,
                            const std::string &datasetName) {
  // Open the HDF5 file in read-only mode
  H5Easy::File file(filename, H5Easy::File::ReadOnly);
  // Load the dataset
  std::vector<int> dataset = H5Easy::load<std::vector<int>>(file, datasetName);
  return dataset;
}

// [[Rcpp::export]]
void replaceH5Dataset(const std::string &filename,
                      const std::string &datasetName,
                      const std::vector<int> &newData) {
  // Open the HDF5 file in read/write mode
  H5Easy::File file(filename, H5Easy::File::ReadWrite);
  // Check if the dataset exists
  if (file.exist(datasetName)) {
    // Delete the existing dataset
    file.unlink(datasetName);
  }
  // Overwrite the dataset with the new data
  H5Easy::dump(file, datasetName, newData, H5Easy::DumpMode::Overwrite);
}
