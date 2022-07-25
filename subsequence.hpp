#ifndef SUBSEQUENCE_H
#define SUBSEQUENCE_H

struct Subsequence {
	double T; // duration
	double C; // accumulated cost
	int W; // delay cost

	// first and last nodes of subsequence
	int first;
	int last;

	inline static Subsequence Concatenate (Subsequence& sigma1, Subsequence& sigma2, double **costMatrix)
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
};

void updateAllSubsequences (Solution *s, double **costMatrix, std::vector<std::vector<Subsequence>>& subseqMatrix);

#endif
