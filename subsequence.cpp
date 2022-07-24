#include "solution.hpp"
#include <vector>

// create header file later

#define SIZE 10

// check later
double t[SIZE][SIZE];

struct Subsequence {
	double T; // duration
	double C; // accumulated cost
	int W; // delay cost

	// first and last nodes of subsequence
	int first;
	int last;

	// why these keywords ?
	inline static Subsequence Concatenate (Subsequence& sigma1, Subsequence& sigma2)
	{
		Subsequence sigma;

		double temp = t[sigma1.last][sigma2.first];

		sigma.W = sigma1.W + sigma2.W;
		sigma.T = sigma1.T + temp + sigma2.T;
		sigma.C = sigma1.C + sigma2.W * (sigma1.T + temp) + sigma2.C;

		sigma.first = sigma1.first;
		sigma.last = sigma2.last;

		return sigma;
	}
};

void updateAllSubsequences (Solution *s, std::vector< std::vector<Subsequence> >& subseq_matrix)
{
	int n = s->sequence.size();


	// subsequences of only one node
	for (int i = 0; i < n; ++i) {
		subseq_matrix[i][i].W = (i > 0);
		subseq_matrix[i][i].C = 0;
		subseq_matrix[i][i].T = 0;
		subseq_matrix[i][i].first = s->sequence[i];
		subseq_matrix[i][i].last = s->sequence[i];
	}

	// remaining subsequences obtained through concatenation
	for (int i = 0; i < n -1; ++i)
	for (int j = i + 1; j < n; ++j) {
		Subsequence currentSubsequence = subseq_matrix[i][j-1];
		Subsequence nodeToAdd = subseq_matrix[j][j];
		subseq_matrix[i][j] = Subsequence::Concatenate(currentSubsequence, nodeToAdd);
	}

	// inverted subsequences (useful for 2-opt)
	for (int i = n -1; i > 0; --i)
	for (int j = i -1; j >= 0; --j) {
		Subsequence currentSubsequence = subseq_matrix[i][j+1];
		Subsequence nodeToAdd = subseq_matrix[j][j];
		subseq_matrix[i][j] = Subsequence::Concatenate(currentSubsequence, nodeToAdd);
	}
}
