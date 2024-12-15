/*
 * Substantial parts of this file were published under MIT license by other
 * authors. Find the original license notice below.
 */

/* MIT License
 *
 * Copyright (c) 2024 sriramlab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "count_data.h"

// [[Rcpp::export]]
int count_samples(std::string filename) {
  ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    i++;
  }
  int n_samples = i - 1;
  return (n_samples);
}

// [[Rcpp::export]]
int count_fam(std::string filename) {
  std::ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    i++;
  }
  return i;
}

// [[Rcpp::export]]
int count_snps_bim(std::string filename) {
  std::ifstream bim_file(filename.c_str(), ios::in);
  if (!bim_file.is_open()) {
    throw std::runtime_error("Error reading .bim file.");
  }
  std::string line;
  int n_snps = 0;
  int linenum = 0;
  while (std::getline(bim_file, line)) {
    linenum++;
    char c = line[0];
    if (c == '#')
      continue;
    if (line.empty())
      continue;
    n_snps++;
  }
  bim_file.close();
  return (n_snps);
}
