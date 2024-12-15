#pragma once
#include <cmath>

struct metaData {
  int wordsize;
  int unitsize;
  unsigned int unitsperword;
  unsigned char mask;
  int ncol;
  int n_samples;
};

metaData set_metadata(int n_samples);
