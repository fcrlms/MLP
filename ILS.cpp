#include <cmath> // INFINITY
#include <vector>

#include "ILS.hpp"
#include "construction.hpp"
#include "localSearch.hpp"
#include "perturb.hpp"
#include "subsequence.hpp"
#include "solution.hpp"

/**
 * Iterated Local Search
 */
Solution ILS (double **matrixAdj, int dimension, int maxIter, int maxIterILS)
{
	Solution bestOfAll;
	bestOfAll.cost = INFINITY;

	auto subseqMatrix = std::vector<std::vector<Subsequence>>(dimension, std::vector<Subsequence>(dimension));

	for (int i = 0; i < maxIter; ++i) {
		Solution s = construction(matrixAdj, dimension);

		updateAllSubsequences(&s, matrixAdj, subseqMatrix);

		s.cost = subseqMatrix[0][dimension-1].C;

		Solution best = s;

		int iterILS = 0;

		while (iterILS <= maxIterILS) {
			localSearch(&s, matrixAdj, subseqMatrix);

			if (s.cost < best.cost) {
				best = s;
				iterILS = 0;
			}

			s = perturb(matrixAdj, &best);

			iterILS += 1;
		}

		if (best.cost < bestOfAll.cost) {
			bestOfAll = best;
		}
	}

	return bestOfAll;
}
