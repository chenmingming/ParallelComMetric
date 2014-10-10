/*
 * Reader.h
 *
 *  Created on: Apr 15, 2014
 *      Author: Miles
 */

#ifndef READER_H_
#define READER_H_

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cassert>
#include <map>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include <vector>

using namespace std;
using namespace std::tr1;

//typedef unordered_map<T, double> innerMap;
//typedef unordered_map<T, innerMap> outerMap;

template<class T>
class Reader {
private:
	string fileName;
public:
	Reader(void);
	Reader(string fileName);

	void setFileName(string fileName);

	string getFileName(void);

	//----------------------For MPIMetric----------------------------------------------//
	// Read community file for each rank
	void getRankCommunity(vector<unordered_set<T> >& rankCommunities, int rank,
			int numprocs);

	// Read communities for each rank and saved as map
	void getRankMapCommunity(unordered_map<T, int>& rankMapCommunities,
			int rank, int numprocs);

	// Read communities for each rank and saved as map and also get their sizes
	void getRankMapCommunity(unordered_map<T, int>& rankMapCommunities,
			unordered_map<int, int>& communitySizes, int rank, int numprocs);

	// Reader detected communities given the real communities and saved as map
	void getDetectedMapCommunity(
			const unordered_map<T, int>& realMapCommunities,
			unordered_map<T, int>& disMapCommunities);

	// Get the network that only contains nodes of the corresponding communities
	double getCommunityNetwork(const unordered_map<T, int>& disMapCommunities,
			unordered_map<T, unordered_map<T, double> >& communityNetwork,
			unordered_set<T>& outCommunityNodes, bool isUnweighted,
			bool isUndirected);

	// Get the network that only contains nodes of the corresponding communities
	double getCommunityReversedNetwork(
			const unordered_map<T, int>& disMapCommunities,
			unordered_map<T, unordered_map<T, double> >& communityInNetwork,
			bool isUnweighted, bool isUndirected);

	void getOutCommunityNodeInfo(unordered_map<int, int>& communitySizes,
			unordered_map<T, int>& disMapCommunities,
			const unordered_set<T>& outCommunityNodes);

	// Get a community from a line of string
	unordered_set<T> getLineCommunity(string lineStr);

	//--------------------------For PthreadMetric---------------------------------------//
	double getNetwork(bool isUnweighted, bool isUndirected,
			unordered_map<T, unordered_map<T, double> >& network);

	double getReversedNetwork(bool isUnweighted, bool isUndirected,
				unordered_map<T, unordered_map<T, double> >& inNet);

	long getCommunity(vector<unordered_set<T> >& communities);

	long getMapCommunity(unordered_map<T, int>& communities);

	~Reader(void);
};

#endif /* READER_H_ */
