/*
 * MPIMetricMain.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Miles
 */

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <mpi.h>
// #include "rdtsc.h"
#include "mm.c"
#include "Reader.cpp"
#include "MPIMetric.cpp"

// #define MAIN_DEBUG
#define GANXIS_MACHINE 1
#define GANXIS_CLOCK_RATE 2100000000
#define GANXIS_CLOCK_RATE_BACKUP 2100000000
#define BLUE_GENE_Q_CLOCK_RATE 1600000000
#define KRATOS_CLOCK_RATE 2667000000

using namespace std;
using namespace std::tr1;

const int OUTPUT_RANK = 0;

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	/* Get the total number of MPI ranks */
	int numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	/* Get the id of this MPI rank */
	int rankid;
	MPI_Comm_rank(MPI_COMM_WORLD, &rankid);

	// The first parameter decide what's kind of metric shall we compute:
	// 1---with ground truth or 0---without ground truth
	int metricType = atoi(argv[1]);
	if (metricType) {
		if (argc != 4) {
			fprintf(stderr,"Usage: %s metricType realCommunityFile detectedCommunityFile\n", argv[0]);
			return 0;
		}
		string realCommunityFile = argv[2];
		string disCommunityFile = argv[3];

#ifdef MAIN_DEBUG
		if (OUTPUT_RANK == rankid) {
			cout << "Real community file = " << realCommunityFile << endl << "Detected community file = "
			<< disCommunityFile << endl;
		}
#endif

		// Reader ground truth community file
		Reader<int> realReader(realCommunityFile);
		// Miles: change set to HashSet
		vector<unordered_set<int> > realCommunities;
		realReader.getRankCommunity(realCommunities, rankid, numprocs);
//		int realComSize = realCommunities.size();
//		cout << "ground truth communities: " << endl;
//		for (int i =0; i < realComSize; ++i) {
//			unordered_set<int> community = realCommunities[i];
//			typename unordered_set<int>::const_iterator comIter;
//			for (comIter = community.begin(); comIter != community.end(); ++comIter) {
//				cout << *comIter << " ";
//			}
//
//			cout << endl;
//		}

		// Read detected community file
		Reader<int> disReader(disCommunityFile);
		// Miles: change set to HashSet
		vector<unordered_set<int> > disCommunities;
		disReader.getRankCommunity(disCommunities, rankid, numprocs);

#ifdef MAIN_DEBUG
		if (OUTPUT_RANK == rankid) {
			int disComSize = disCommunities.size();
			cout << disComSize << " discovered communities: " << endl;
			for (int i = 0; i < disComSize; ++i) {
				unordered_set<int> community = disCommunities[i];
				typename unordered_set<int>::const_iterator comIter;
				for (comIter = community.begin(); comIter != community.end(); ++comIter) {
					cout << *comIter << " ";
				}
				cout << endl;
			}
		}
#endif

		MPIMetric<int> mpiMetric;

		// Compute Information Entropy metrics
		double entropyMetrics[2];
		double entropyCycles[2] = {0.0};
		double startInfoEntropyCycle = rdtsc();
		mpiMetric.computeInfoEntropyMetric(rankid, numprocs, realCommunities, disCommunities, entropyMetrics, entropyCycles);
		double endInfoEntropyCycle = rdtsc();
		double rankEntropyExeCyle = endInfoEntropyCycle - startInfoEntropyCycle;
		double avgEntropyExeCycle, avgEntropyCompCycle, avgEntropyMsgCycle;
		MPI_Allreduce(&rankEntropyExeCyle, &avgEntropyExeCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&entropyCycles[0], &avgEntropyCompCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&entropyCycles[1], &avgEntropyMsgCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		avgEntropyExeCycle /= numprocs;
		avgEntropyCompCycle /= numprocs;
		avgEntropyMsgCycle /= numprocs;

		double entropyExeTime, entropyCompTime, entropyMsgTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			entropyExeTime = avgEntropyExeCycle / GANXIS_CLOCK_RATE;
			entropyCompTime = avgEntropyCompCycle / GANXIS_CLOCK_RATE;
			entropyMsgTime = avgEntropyMsgCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			entropyExeTime = avgEntropyExeCycle / BLUE_GENE_Q_CLOCK_RATE;
			entropyCompTime = avgEntropyCompCycle / BLUE_GENE_Q_CLOCK_RATE;
			entropyMsgTime = avgEntropyMsgCycle / BLUE_GENE_Q_CLOCK_RATE;
		}

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << numprocs << "\t" << entropyExeTime << "\t" << entropyCompTime << "\t" << entropyMsgTime << endl;
		}

		if (rankid == OUTPUT_RANK) {
			cout << "MAIN: VI = " << entropyMetrics[0] << ", NMI = " << entropyMetrics[1] << endl;
		}
#endif

		// Compute Clustering Matching metrics
		double clusteringMetrics[2];
		double clusteringCycles[2] = {0.0};
		double startClusteringCycle = rdtsc();
		mpiMetric.computeClusterMatchingMetric(rankid, numprocs, realCommunities, disCommunities, clusteringMetrics, clusteringCycles);
		double endClusteringCycle = rdtsc();
		double rankClusteringExeCyle = endClusteringCycle - startClusteringCycle;
		double avgClusteringExeCycle, avgClusteringCompCycle, avgClusteringMsgCycle;
		MPI_Allreduce(&rankClusteringExeCyle, &avgClusteringExeCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&clusteringCycles[0], &avgClusteringCompCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&clusteringCycles[1], &avgClusteringMsgCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		avgClusteringExeCycle /= numprocs;
		avgClusteringCompCycle /= numprocs;
		avgClusteringMsgCycle /= numprocs;

		double clusteringExeTime, clusteringCompTime, clusteringMsgTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			clusteringExeTime = avgClusteringExeCycle / GANXIS_CLOCK_RATE;
			clusteringCompTime = avgClusteringCompCycle / GANXIS_CLOCK_RATE;
			clusteringMsgTime = avgClusteringMsgCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			clusteringExeTime = avgClusteringExeCycle / BLUE_GENE_Q_CLOCK_RATE;
			clusteringCompTime = avgClusteringCompCycle / BLUE_GENE_Q_CLOCK_RATE;
			clusteringMsgTime = avgClusteringMsgCycle / BLUE_GENE_Q_CLOCK_RATE;
		}

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << numprocs << "\t" << clusteringExeTime << "\t" << clusteringCompTime << "\t" << clusteringMsgTime << endl;
		}

		if (rankid == OUTPUT_RANK) {
			cout << "MAIN: FMeasure = " << clusteringMetrics[0] << ", NVD = " << clusteringMetrics[1] << endl;
		}
#endif

		// Compute Index metrics
		unordered_map<int, int> realMapCommunities;
		realReader.getRankMapCommunity(realMapCommunities, rankid, numprocs);

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << realMapCommunities.size() << " real node communities:" << endl;
			typename unordered_map<int, int>::const_iterator realIter;
			for (realIter = realMapCommunities.begin(); realIter != realMapCommunities.end(); ++realIter) {
				cout << realIter->first << "\t" << realIter->second << endl;
			}
		}
#endif

		unordered_map<int, int> disMapCommunities;
		disReader.getDetectedMapCommunity(realMapCommunities, disMapCommunities);

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << disMapCommunities.size() << " detected node communities:" << endl;
			typename unordered_map<int, int>::const_iterator disIter;
			for (disIter = disMapCommunities.begin(); disIter != disMapCommunities.end(); ++disIter) {
				cout << disIter->first << "\t" << disIter->second << endl;
			}
		}
#endif
		double indexMetrics[3];
		double indexCycles[2] = {0.0};
		double startIndexCycle = rdtsc();
		mpiMetric.computeIndexMetric(rankid, numprocs,realMapCommunities, disMapCommunities,indexMetrics, indexCycles);
		double endIndexCycle = rdtsc();
		double rankIndexExeCyle = endIndexCycle - startIndexCycle;
		double avgIndexExeCycle, avgIndexCompCycle, avgIndexMsgCycle;
		MPI_Allreduce(&rankIndexExeCyle, &avgIndexExeCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&indexCycles[0], &avgIndexCompCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(&indexCycles[1], &avgIndexMsgCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		avgIndexExeCycle /= numprocs;
		avgIndexCompCycle /= numprocs;
		avgIndexMsgCycle /= numprocs;

		double indexExeTime, indexCompTime, indexMsgTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			indexExeTime = avgIndexExeCycle / GANXIS_CLOCK_RATE;
			indexCompTime = avgIndexCompCycle / GANXIS_CLOCK_RATE;
			indexMsgTime = avgIndexMsgCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			indexExeTime = avgIndexExeCycle / BLUE_GENE_Q_CLOCK_RATE;
			indexCompTime = avgIndexCompCycle / BLUE_GENE_Q_CLOCK_RATE;
			indexMsgTime = avgIndexMsgCycle / BLUE_GENE_Q_CLOCK_RATE;
		}

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << numprocs << "\t" << indexExeTime << "\t" << indexCompTime << "\t" << indexMsgTime << endl;
		}

		if (rankid == OUTPUT_RANK) {
			cout << "MAIN: RI = " << indexMetrics[0] << ", ARI = " << indexMetrics[1] << ", JI = " << indexMetrics[2] << endl;
		}
#endif

		if(rankid == OUTPUT_RANK) {
			cout << numprocs << "\t" << entropyExeTime << "\t" << entropyCompTime << "\t" << entropyMsgTime << "\t"
			<< clusteringExeTime << "\t" << clusteringCompTime << "\t" << clusteringMsgTime << "\t"
			<< indexExeTime << "\t" << indexCompTime << "\t" << indexMsgTime << endl;

			cout << numprocs << "\t" << entropyMetrics[0] << "\t" << entropyMetrics[1] << "\t"
			<< clusteringMetrics[0] << "\t" << clusteringMetrics[1] << "\t" << indexMetrics[0] << "\t"
			<< indexMetrics[1] << "\t" << indexMetrics[2] << endl;
		}
	}
	// For metric without ground truth communities
	else {
		if (argc < 4) {
			fprintf(stderr,"Usage: %s metricType detectedCommunityFile networkFile [isUnweighted isUndirected]\n", argv[0]);
			return 0;
		}
		string disCommunityFile = argv[2];
		string networkFile = argv[3];
		bool isUnweighted = true;
		if (argc > 4) {
			isUnweighted = atoi(argv[4]);
		}
		bool isUndirected = true;
		if (argc > 5) {
			isUndirected = atoi(argv[5]);
		}

#ifdef MAIN_DEBUG
		if (OUTPUT_RANK == rankid) {
			cout << "Detected community file = " << disCommunityFile << endl << "Network file = " << networkFile << endl
			<< "isUnweighted = " << isUnweighted << ", isUndirected = " << isUndirected << endl;
		}
#endif

		// Read detected community file
		Reader<int> disReader(disCommunityFile);
		// Miles: change set to HashSet
		unordered_map<int, int> communitySizes;
		unordered_map<int, int> disMapCommunities;
		disReader.getRankMapCommunity(disMapCommunities, communitySizes, rankid, numprocs);

		// Read the network that only contains nodes of the corresponding communities
		unordered_map<int, unordered_map<int, double> > communityNetwork;
		unordered_set<int> outCommunityNodes;
		Reader<int> netReader(networkFile);
		double comNetWeight = netReader.getCommunityNetwork(disMapCommunities, communityNetwork, outCommunityNodes, isUnweighted, isUndirected);

#ifdef MAIN_DEBUG
		if (OUTPUT_RANK == rankid) {
			cout << "community node size = " << disMapCommunities.size() << ", net node size = " << communityNetwork.size() << endl;
			typename unordered_map<int, unordered_map<int, double> >::const_iterator netIter;
			for (netIter = communityNetwork.begin(); netIter != communityNetwork.end(); ++netIter) {
				int srcId = netIter->first;
				unordered_map<int, double> nbs = netIter->second;

				cout << "MAIN: nodeId = " << srcId << ", nbSize = " << nbs.size() << endl;
			}
		}
#endif

		disReader.getOutCommunityNodeInfo(communitySizes, disMapCommunities, outCommunityNodes);
		outCommunityNodes.clear();

		MPIMetric<int> mpiMetric;
		double metrics[8];
		double startExeCycle = rdtsc();
		mpiMetric.computeMetricWithoutGroundTruth(rankid, comNetWeight, communitySizes, disMapCommunities, communityNetwork, isUndirected, metrics);
		double endExeCycle = rdtsc();
		double rankExeCycle = endExeCycle - startExeCycle;
		double avgExeCycle;
		MPI_Allreduce(&rankExeCycle, &avgExeCycle, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		avgExeCycle /= numprocs;
		double exeTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			exeTime = avgExeCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			exeTime = avgExeCycle / BLUE_GENE_Q_CLOCK_RATE;
		}

#ifdef MAIN_DEBUG
		if (rankid == OUTPUT_RANK) {
			cout << "MAIN: Q = " << metrics[0] << ", Qds = " << metrics[1] << ", intraEdges = " << metrics[2]
			<< ", intraDensity = " << metrics[3] << ", contraction = " << metrics[4] << ", interEdges = " << metrics[5]
			<< ", expansion = " << metrics[6] << ", conductance = " << metrics[7] << endl;
		}
#endif

		if (rankid == OUTPUT_RANK) {
			cout << numprocs << "\t" << exeTime << endl;

			cout << numprocs << "\t" << metrics[0] << "\t" << metrics[1] << "\t" << metrics[2] << "\t"
			<< metrics[3] << "\t" << metrics[4] << "\t" << metrics[5] << "\t"
			<< metrics[6] << "\t" << metrics[7] << endl;
		}
	}

	MPI_Finalize();
	return 0;
}
