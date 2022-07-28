#include <cstdlib> // rand()
#include <vector> // swap()
#include <cmath> // ceil()

#include "localSearch.hpp"
#include "subsequence.hpp"
#include "solution.hpp"

#define NL_SIZE 5

bool bestImprovementSwap (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix);
bool bestImprovement2Opt (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix);
bool bestImprovementOrOpt (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix, int n);

/**
 * Based on RVND (Random Variable Neighborhood Descent)
 */
void localSearch (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix)
{
	int NL[NL_SIZE] = { 1, 2, 3, 4, 5 };
	int offset = 0;

	while (offset < NL_SIZE -1) {
		bool improved = false;

		int randomIndex = std::rand() % (NL_SIZE - offset);
		int n = NL[randomIndex];

		switch (n) {
		case 1:
			improved = bestImprovementSwap(s, matrixAdj, subseqMatrix);
			break;
		case 2:
			improved = bestImprovement2Opt(s, matrixAdj, subseqMatrix);
			break;
		case 3:
			improved = bestImprovementOrOpt(s, matrixAdj, subseqMatrix, 1); // Reinsertion
			break;
		case 4:
			improved = bestImprovementOrOpt(s, matrixAdj, subseqMatrix, 2); // Or-opt2
			break;
		case 5:
			improved = bestImprovementOrOpt(s, matrixAdj, subseqMatrix, 3); // Or-opt3
			break;
		}

		if (improved) {
			// reset
			offset = 0;
		} else {
			// remove choice
			std::swap(NL[randomIndex], NL[NL_SIZE -1 -offset]);
			offset += 1;
		}
	}

	return;
}

/**
 * @param n number of nodes
 * @param m cost matrix
 * @param sM subsequence matrix
 * @param i ith node position
 * @param j jth node position
 */
double calculateSwapCost (int n, double **m, std::vector<std::vector<Subsequence>>& sM, int i, int j)
{
	Subsequence old = sM[0][n-1];

	Subsequence sigma = sM[0][i-1];

	sigma = Concatenate(sigma, sM[j][j], m);

	if (i + 1 != j) {
		sigma = Concatenate(sigma, sM[i+1][j-1], m);
	}

	sigma = Concatenate(sigma, sM[i][i], m);

	sigma = Concatenate(sigma, sM[j+1][n-1], m);

	return sigma.C - old.C;
}

bool bestImprovementSwap (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix)
{
	double bestDelta = 0;
	double best_i = 0;
	double best_j = 0;

	int range = s->sequence.size() -1;

	for (int i = 1; i < range -1; ++i)
	for (int j = i +1; j < range; ++j) {
		double thisDelta = calculateSwapCost(range +1, matrixAdj, subseqMatrix, i, j);

		if (thisDelta < bestDelta) {
			bestDelta = thisDelta;
			best_i = i;
			best_j = j;
		}
	}

	if (bestDelta < 0) {
		std::swap(s->sequence[best_i], s->sequence[best_j]);

		s->cost += bestDelta;

		updateAllSubsequences(s, matrixAdj, subseqMatrix);

		return true;
	}

	return false;
}

/**
 * @param n number of nodes
 * @param m cost matrix
 * @param sM subsequence matrix
 * @param i subsequence start
 * @param j subsequence end
 */
double calculate2OptCost (int n, double **m, std::vector<std::vector<Subsequence>>& sM, int i, int j)
{
	Subsequence old = sM[0][n-1];

	Subsequence sigma = sM[0][i];

	// the subsequence from i+1 to j is inverted and reinserted
	sigma = Concatenate(sigma, sM[j][i+1], m);

	sigma = Concatenate(sigma, sM[j+1][n-1], m);

	return sigma.C - old.C;
}

void exec2Opt (std::vector<int>& s, int i, int j)
{
	// node i+1 to node j
	int range = std::ceil((j - i) / 2);

	// inverts the segment between the edges {i, i+1} and {j, j+1}
	for (int n = 0; n < range; ++n) {
		int node_i = i+1 + n; // i+1 is the initial node
		int node_j = j -n; // j is the final node

		std::swap(s[node_i], s[node_j]);
	}
}

bool bestImprovement2Opt (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix)
{
	double bestDelta = 0;

	double best_i = 0; // ith edge
	double best_j = 0; // jth edge

	int range = s->sequence.size() -1;

	/**
	 * When i=0 we have a special case where the final edge can't be chosen
	 * because it is adjacent to the first edge since it is a loop
	 */
	for (int j = 2; j < range; ++j) {
		double thisDelta = calculate2OptCost(range +1, matrixAdj, subseqMatrix, 0, j);

		if (thisDelta < bestDelta) {
			bestDelta = thisDelta;
			best_i = 0;
			best_j = j;
		}
	}

	// remaining cases
	for (int i = 1; i < range -2; ++i)
	for (int j = i +2; j < range; ++j) {
		double thisDelta = calculate2OptCost(range +1, matrixAdj, subseqMatrix, i, j);

		if (thisDelta < bestDelta) {
			bestDelta = thisDelta;
			best_i = i;
			best_j = j;
		}
	}

	if (bestDelta < 0) {
		exec2Opt(s->sequence, best_i, best_j);

		s->cost += bestDelta;

		updateAllSubsequences(s, matrixAdj, subseqMatrix);

		return true;
	}

	return false;
}

void execOrOpt(std::vector<int>& s, int best_start, int n, int best_posInsert)
{
	// sequence that will be repositioned
	int chain[3];

	for (int i = 0; i < n; ++i) {
		chain[i] = s[best_start + i];
	}

	bool forward = (best_posInsert >= best_start);

	if (forward) {
		// moving elements ahead of the chain n positions backwards
		for (int i = best_start + n; i <= best_posInsert; ++i) {
			s[i - n] = s[i];
		}

		/**
		 * The insertion point has moved n positions backwards
		 * so we begin at (best_posInsert - n + 1)
		 */
		for (int i = 0; i < n; ++i) {
			s[(best_posInsert - n +1) + i] = chain[i];
		}
	} else {
		/**
		 * Moving elements behind the chain (until the insertion point)
		 * n positions forward
		 */
		for (int i = best_start -1; i > best_posInsert; --i) {
			s[i + n] = s[i];
		}

		/**
		 * The insertion point hasn't moved so we begin inserting at
		 * (best_posInsert +1)
		 */
		for (int i = 0; i < n; ++i) {
			s[(best_posInsert +1) + i] = chain[i];
		}
	}
}

/**
 * @param n number of nodes
 * @param m cost matrix
 * @param sM subsequence matrix
 * @param start start of chain
 * @param end end of chain (last element of chain is end-1)
 * @param pos position to insert the chain (first node of chain will be at pos+1)
 */
double calculateOrOptCost (int n, double **m, std::vector<std::vector<Subsequence>>& sM, int start, int end, int pos)
{
	Subsequence old = sM[0][n-1];

	// position to insert comes before the chain
	if (pos < start) {
		Subsequence sigma = sM[0][pos];
		sigma = Concatenate(sigma, sM[start][end-1], m);
		sigma = Concatenate(sigma, sM[pos+1][start-1], m);
		sigma = Concatenate(sigma, sM[end][n-1], m);
		return sigma.C - old.C;
	}
	// position to insert comes after the chain
	else {
		Subsequence sigma = sM[0][start-1];
		sigma = Concatenate(sigma, sM[end][pos], m);
		sigma = Concatenate(sigma, sM[start][end-1], m);
		sigma = Concatenate(sigma, sM[pos+1][n-1], m);
		return sigma.C - old.C;
	}
}

/**
 * @param n size of chain
 */
bool bestImprovementOrOpt (Solution *s, double **matrixAdj, std::vector<std::vector<Subsequence>>& subseqMatrix, int n)
{
	double bestDelta = 0;
	double best_start = 0;
	double best_posInsert = 0;

	int range = s->sequence.size() -n;
	int totalNodes = s->sequence.size();

	/**
	 * The first element of the chain is start
	 * The last element of the chain is end-1
	 */
	for (int start = 1; start < range; ++start) {
		int end = start + n;

		// i is the position to insert, the first node will be put at i+1
		for (int i = 0; i < start -1; ++i) {
			double thisDelta = calculateOrOptCost(totalNodes, matrixAdj, subseqMatrix, start, end, i);

			if (thisDelta < bestDelta) {
				bestDelta = thisDelta;
				best_start = start;
				best_posInsert = i;
			}
		}

		for (int i = end; i < range; ++i) {
			double thisDelta = calculateOrOptCost(totalNodes, matrixAdj, subseqMatrix, start, end, i);

			if (thisDelta < bestDelta) {
				bestDelta = thisDelta;
				best_start = start;
				best_posInsert = i;
			}
		}
	}

	if (bestDelta < 0) {
		execOrOpt(s->sequence, best_start, n, best_posInsert);

		s->cost += bestDelta;

		updateAllSubsequences(s, matrixAdj, subseqMatrix);

		return true;
	}

	return false;
}
