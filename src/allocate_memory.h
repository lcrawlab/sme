//
// Created by Julian Stamp on 3/29/24.
//

#ifndef MMER_ALLOCATE_MEMORY_H
#define MMER_ALLOCATE_MEMORY_H

#include "genotype.h"
#include "mme.h"

void allocate_memory(int n_randvecs, const genotype &genotype_block,
                     double *&partialsums, double *&sum_op, double *&yint_e,
                     double *&yint_m, double **&y_e, double **&y_m);

void deallocate_memory(double *partialsums, double *sum_op, double *yint_e,
                       double *yint_m, double **y_e, double **y_m,
                       const genotype &genotype_block);

#endif // MMER_ALLOCATE_MEMORY_H
