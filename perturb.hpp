#ifndef PERTURB_H
#define PERTURB_H

#include "subsequence.hpp"
#include "solution.hpp"

Solution perturb (Solution *best, double **costMatrix, std::vector<std::vector<Subsequence>>& subseqMatrix);

#endif
