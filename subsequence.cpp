#include <vector>

#include "subsequence.hpp"
#include "solution.hpp"

// later create a specific function to update the subsequences according to what movement was made
void updateAllSubsequences (Solution *s, double **costMatrix, std::vector<std::vector<Subsequence>>& subseqMatrix)
{
	int n = s->sequence.size();

	// subsequences of only one node
	for (int i = 0; i < n; ++i) {
		subseqMatrix[i][i].W = (i > 0);
		subseqMatrix[i][i].C = 0;
		subseqMatrix[i][i].T = 0;
		subseqMatrix[i][i].first = s->sequence[i];
		subseqMatrix[i][i].last = s->sequence[i];
	}

	// remaining subsequences obtained through concatenation
	for (int i = 0; i < n -1; ++i)
	for (int j = i +1; j < n; ++j) {
		Subsequence currentSubsequence = subseqMatrix[i][j-1];
		Subsequence nodeToAdd = subseqMatrix[j][j];
		subseqMatrix[i][j] = Concatenate(currentSubsequence, nodeToAdd, costMatrix);
	}

	// inverted subsequences (useful for 2-opt)
	for (int i = n -1; i > 0; --i)
	for (int j = i -1; j >= 0; --j) {
		Subsequence currentSubsequence = subseqMatrix[i][j+1];
		Subsequence nodeToAdd = subseqMatrix[j][j];
		subseqMatrix[i][j] = Concatenate(currentSubsequence, nodeToAdd, costMatrix);
	}
}

Subsequence Concatenate (Subsequence& sigma1, Subsequence& sigma2, double **costMatrix)
{
	Subsequence sigma;

	double temp = costMatrix[sigma1.last][sigma2.first];

	sigma.W = sigma1.W + sigma2.W;
	sigma.T = sigma1.T + temp + sigma2.T;
	sigma.C = sigma1.C + sigma2.W * (sigma1.T + temp) + sigma2.C;

	sigma.first = sigma1.first;
	sigma.last = sigma2.last;

	return sigma;
}
