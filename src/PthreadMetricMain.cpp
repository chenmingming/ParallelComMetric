/*
 * PthreadMetricMain.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: Miles
 */

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <pthread.h>
// #include "rdtsc.h"
#include "mm.c"
#include "PthreadMetric.cpp"

// #define MAIN_DEBUG
#define GANXIS_MACHINE 1
#define GANXIS_CLOCK_RATE 2100000000
#define GANXIS_CLOCK_RATE_BACKUP 2100000000
#define BLUE_GENE_Q_CLOCK_RATE 1600000000
#define KRATOS_CLOCK_RATE 2667000000

using namespace std;
using namespace std::tr1;

int main(int argc, char **argv) {
	// The first parameter decide what's kind of metric shall we compute:
	// 1---with ground truth or 0---without ground truth
	int numThreads = atoi(argv[1]);
	int metricType = atoi(argv[2]);
	if (metricType) {
		if (argc != 5) {
			fprintf(stderr,"Usage: %s numThreads metricType realCommunityFile detectedCommunityFile\n", argv[0]);
			return 0;
		}
		string realCommunityFile = argv[3];
		string disCommunityFile = argv[4];

		PthreadMetric<int> pthreadMetric(numThreads, realCommunityFile,
				disCommunityFile);
		double entropyExeCycle = pthreadMetric.computeInfoEntropyMetric();
		double entropyExeTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			entropyExeTime = entropyExeCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			entropyExeTime = entropyExeCycle / BLUE_GENE_Q_CLOCK_RATE;
		}
		double VI = pthreadMetric.getVI();
		double NMI = pthreadMetric.getNMI();

#ifdef MAIN_DEBUG
		cout << "MAIN: VI = " << VI << ", NMI = " << NMI << endl;
#endif

		double clusteringExeCycle = pthreadMetric.computeClusterMatchingMetric();
		double clusteringExeTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			clusteringExeTime = clusteringExeCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			clusteringExeTime = clusteringExeCycle / BLUE_GENE_Q_CLOCK_RATE;
		}
		double fMeasure = pthreadMetric.getFMeasure();
		double NVD = pthreadMetric.getNVD();

#ifdef MAIN_DEBUG
		cout << "MAIN: fMeasure = " << fMeasure << ", NVD = " << NVD << endl;
#endif

		double indexExeCycle = pthreadMetric.computeIndexMetric();
		double indexExeTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			indexExeTime = indexExeCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			indexExeTime = indexExeCycle / BLUE_GENE_Q_CLOCK_RATE;
		}
		double RI = pthreadMetric.getRI();
		double ARI = pthreadMetric.getARI();
		double JI = pthreadMetric.getJI();

#ifdef MAIN_DEBUG
		cout << "MAIN: RI = " << RI << ", ARI = " << ARI << ", JI = " << JI << endl;
#endif

		cout << numThreads << "\t" << entropyExeTime << "\t" << clusteringExeTime << "\t" << indexExeTime << endl;
		cout << numThreads << "\t" << VI << "\t" << NMI << "\t" << fMeasure << "\t" << NVD << "\t" << RI << "\t" << ARI << "\t" << JI << endl;
	}
	// For metric without ground truth communities
	else {
		if (argc < 5) {
			fprintf(stderr,"Usage: %s numThreads metricType detectedCommunityFile networkFile [isUnweighted isUndirected]\n", argv[0]);
			return 0;
		}
		string disCommunityFile = argv[3];
		string networkFile = argv[4];
		bool isUnweighted = true;
		if (argc > 5) {
			isUnweighted = atoi(argv[5]);
		}
		bool isUndirected = true;
		if (argc > 6) {
			isUndirected = atoi(argv[6]);
		}

		PthreadMetric<int> pthreadMetric(numThreads, networkFile, isUnweighted,
				isUndirected, disCommunityFile);
		double exeCycle = pthreadMetric.computeMetricWithoutGroundTruth();
		double exeTime;
		if (GANXIS_MACHINE) {
			// ganxis machine
			exeTime = exeCycle / GANXIS_CLOCK_RATE;
		} else {
			// Blue Gene/Q
			exeTime = exeCycle / BLUE_GENE_Q_CLOCK_RATE;
		}
		double modularity = pthreadMetric.getModularity();
		double Qds = pthreadMetric.getQds();
		double intraEdges = pthreadMetric.getIntraEdges();
		double intraDensity = pthreadMetric.getIntraDensity();
		double contraction = pthreadMetric.getContraction();
		double interEdges = pthreadMetric.getInterEdges();
		double expansion = pthreadMetric.getExpansion();
		double conductance = pthreadMetric.getConductance();

#ifdef MAIN_DEBUG
		cout << "MAIN: Q = " << modularity << ", Qds = " << Qds << ", intraEdges = " << intraEdges
		<< ", intraDensity = " << intraDensity << ", contraction = " << contraction << ", interEdges = " << interEdges
		<< ", expansion = " <<expansion << ", conductance = " << conductance << endl;
#endif

		cout << numThreads << "\t" << exeTime << endl;
		cout << numThreads << "\t" << modularity << "\t" << Qds << "\t" << intraEdges << "\t"
		<< intraDensity << "\t" << contraction << "\t" << interEdges << "\t" << expansion
		<< "\t" << conductance << endl;
	}

	return 0;
}

