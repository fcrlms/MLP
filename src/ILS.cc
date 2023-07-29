#include <cmath> // INFINITY
#include <vector>

#include "ILS.h"
#include "construction.h"
#include "localSearch.h"
#include "perturb.h"
#include "subsequence.h"
#include "solution.h"

/**
 * Iterated Local Search
 */
Solution ILS (double **matrixAdj, int dimension, int maxIter, int maxIterILS)
{
	Solution bestOfAll;
	bestOfAll.cost = INFINITY;

	// the first node repeats as the last node and every other node appears only once
	// then (dimension +1) is needed to represent all subsequences
	auto subseqMatrix = std::vector<std::vector<Subsequence>>(dimension +1, std::vector<Subsequence>(dimension +1));

	for (int i = 0; i < maxIter; ++i) {
		Solution s = construction(matrixAdj, dimension);

		updateAllSubsequences(&s, matrixAdj, subseqMatrix);

		s.cost = subseqMatrix[0][dimension].C;

		Solution best = s;

		int iterILS = 0;

		while (iterILS <= maxIterILS) {
			localSearch(&s, matrixAdj, subseqMatrix);

			if (s.cost < best.cost) {
				best = s;
				iterILS = 0;
			}

			s = perturb(&best, matrixAdj, subseqMatrix);

			iterILS += 1;
		}

		if (best.cost < bestOfAll.cost) {
			bestOfAll = best;
		}
	}

	return bestOfAll;
}
