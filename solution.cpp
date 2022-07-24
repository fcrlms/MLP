#include <iostream>
#include "solution.hpp"

void printSolution (Solution *s) {
	int sequenceSize = s->sequence.size();

	for (int i = 0; i < sequenceSize -1; ++i) {
		std::cout << s->sequence[i] << " -> ";
	}

	std::cout << s->sequence.back() << std::endl;
}

void updateSolutionCost (Solution *s, double **costMatrix) {
	double cost = 0;
	int sequenceSize = s->sequence.size();

	for (int i = 0; i < sequenceSize -1; ++i) {
		int currNode = s->sequence[i];
		int nextNode = s->sequence[i +1];

		/**
		 * The accumulated cost is a summation and the ith cost (costMatrix[currNode][nextNode])
		 * will repeat (sequenceSize - 1 - i) times
		 */
		cost += (sequenceSize - 1 - i) * costMatrix[currNode][nextNode];
	}

	s->cost = cost;
}
