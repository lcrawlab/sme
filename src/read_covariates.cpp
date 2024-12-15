#include "read_covariates.h"

int read_covariates(const bool &standardize, const int &Nind,
                    const std::string &filename, MatrixXdr &covariate) {
  std::string covname = "";
  ifstream ifs(filename.c_str(), ios::in);
  std::string line;
  std::istringstream in;
  int covIndex = 0;
  std::getline(ifs, line);
  in.str(line);
  string b;
  vector<vector<int>> missing;
  int covNum = 0;
  while (in >> b) {
    if (b != "FID" && b != "IID") {
      missing.push_back(vector<int>()); // push an empty row
      if (b == covname && covname != "")
        covIndex = covNum;
      covNum++;
    }
  }
  vector<double> cov_sum(covNum, 0);
  if (covname == "") {
    covariate.resize(Nind, covNum);
  } else {
    covariate.resize(Nind, 1);
  }

  int j = 0;
  while (std::getline(ifs, line)) {
    in.clear();
    in.str(line);
    string temp;
    in >> temp;
    in >> temp; // FID IID
    for (int k = 0; k < covNum; k++) {

      in >> temp;
      if (temp == "NA") {
        missing[k].push_back(j);
        continue;
      }
      double cur = atof(temp.c_str());
      if (cur == -9) {
        missing[k].push_back(j);
        continue;
      }
      if (covname == "") {
        cov_sum[k] = cov_sum[k] + cur;
        covariate(j, k) = cur;
      } else if (k == covIndex) {
        covariate(j, 0) = cur;
        cov_sum[k] = cov_sum[k] + cur;
      }
    }
    j++;
  }
  // compute cov mean and impute
  for (int a = 0; a < covNum; a++) {
    int missing_num = missing[a].size();
    cov_sum[a] = cov_sum[a] / (Nind - missing_num);

    for (int b = 0; b < missing_num; b++) {
      int index = missing[a][b];
      if (covname == "")
        covariate(index, a) = cov_sum[a];
      else if (a == covIndex)
        covariate(index, 0) = cov_sum[a];
    }
  }
  if (standardize) {
    MatrixXdr cov_std;
    cov_std.resize(1, covNum);
    MatrixXdr sum = covariate.colwise().sum();
    MatrixXdr sum2 = (covariate.cwiseProduct(covariate)).colwise().sum();
    MatrixXdr temp;
    for (int b = 0; b < covNum; b++) {
      cov_std(0, b) = sum2(0, b) + Nind * cov_sum[b] * cov_sum[b] -
                      2 * cov_sum[b] * sum(0, b);
      cov_std(0, b) = sqrt((Nind - 1) / cov_std(0, b));
      double scalar = cov_std(0, b);
      for (int j = 0; j < Nind; j++) {
        covariate(j, b) = covariate(j, b) - cov_sum[b];
        covariate(j, b) = covariate(j, b) * scalar;
      }
    }
  }
  return covNum;
}
