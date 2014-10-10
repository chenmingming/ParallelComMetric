/*
 * MPIMetric.h
 *
 *  Created on: Apr 16, 2014
 *      Author: Miles
 */

#ifndef MPIMETRIC_H_
#define MPIMETRIC_H_

#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <tr1/unordered_set>
// #include "rdtsc.h"
#include "mm.c"
#include <mpi.h>

// #define MPI_DEBUG

using namespace std;
using namespace std::tr1;

template<class T>
class MPIMetric {
private:
	static const int END_OF_COMMUNITY = -1;
	static const int OUTPUT_RANK = 0;

	int getIntersectionNum(unordered_set<T>& scomm, unordered_set<T>& gcomm);
	T** malloc2d(int row, int column);
	void free2d(T **arr);
public:
	MPIMetric(void);

	void computeMetricWithoutGroundTruth(int rankid, double comNetWeight,
			unordered_map<int, int>& communitySizes,
			unordered_map<T, int>& disMapCommunities,
			const unordered_map<T, unordered_map<T, double> >& communityNetwork,
			const unordered_map<T, unordered_map<T, double> >& communityInNetwork,
			bool isUndirected, double *metrics);

	void computeInfoEntropyMetric(int rankid, int numprocs,
			const vector<unordered_set<T> >& realCommunities,
			const vector<unordered_set<T> >& disCommunities, double *metrics,
			double *entropyCycles);

	void computeInfoEntropyMetric_backup(int rankid, int numprocs,
			const vector<unordered_set<T> >& realCommunities,
			const vector<unordered_set<T> >& disCommunities, double *metrics);

	void computeClusterMatchingMetric(int rankid, int numprocs,
			const vector<unordered_set<T> >& realCommunities,
			const vector<unordered_set<T> >& disCommunities, double *metrics,
			double *clusteringCycles);

	void computeClusterMatchingMetric_backup(int rankid, int numprocs,
			const vector<unordered_set<T> >& realCommunities,
			const vector<unordered_set<T> >& disCommunities, double *metrics);

	void computeIndexMetric(int rankid, int numprocs,
			unordered_map<T, T>& realMapCommunities,
			unordered_map<T, T>& disMapCommunities, double *metrics,
			double *indexCycles);

	void computeIndexMetric_backup(int rankid, int numprocs,
			unordered_map<T, T>& realMapCommunities,
			unordered_map<T, T>& disMapCommunities, double *metrics,
			double *indexCycles);

	~MPIMetric(void);
};

#endif /* MPIMETRIC_H_ */
