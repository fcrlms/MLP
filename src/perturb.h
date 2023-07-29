#ifndef PERTURB_H
#define PERTURB_H

#include "subsequence.h"
#include "solution.h"

Solution perturb (Solution *best, double **costMatrix, std::vector<std::vector<Subsequence>>& subseqMatrix);

#endif
