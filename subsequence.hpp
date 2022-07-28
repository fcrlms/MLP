#ifndef SUBSEQUENCE_H
#define SUBSEQUENCE_H

#include "solution.hpp"

struct Subsequence {
	double T; // duration
	double C; // accumulated cost
	int W; // delay cost

	// first and last nodes of subsequence
	int first;
	int last;
};

void updateAllSubsequences (Solution *s, double **costMatrix, std::vector<std::vector<Subsequence>>& subseqMatrix);
Subsequence Concatenate (Subsequence& sigma1, Subsequence& sigma2, double **costMatrix);
void Append (Subsequence* sigma1, Subsequence& sigma2, double **costMatrix);

#endif
