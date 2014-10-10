/*
 * Reader.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Miles
 */

#ifndef READER_CPP
#define READER_CPP

#include "Reader.h"

template<class T>
Reader<T>::Reader(void) {
	this->fileName = "";
}

template<class T>
Reader<T>::Reader(string fileName) {
	this->fileName = fileName;
}

template<class T>
void Reader<T>::setFileName(string fileName) {
	this->fileName = fileName;
}

template<class T>
string Reader<T>::getFileName(void) {
	return this->fileName;
}

template<class T>
void Reader<T>::getRankCommunity(vector<unordered_set<T> >& rankCommunities,
		int rank, int numprocs) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getRankCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return;
	}

	// Read until all the lines are read
	int count = 0;
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		// Each rank read its own communities. This method could also achieve rough loading balance.
		if (count % numprocs == rank) {
			unordered_set<T> community = this->getLineCommunity(lineStr);
//			if (community.empty()) {
//				continue;
//			}
			rankCommunities.push_back(community);
		}

		++count;
	}

	communityFp.close();
}

template<class T>
void Reader<T>::getRankMapCommunity(unordered_map<T, int>& rankMapCommunities,
		int rank, int numprocs) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getRankMapCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return;
	}

	// Community id starts from 0
	int count = 0;
	// Read until all the lines are read
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		// Each rank read its own communities. This method could also achieve rough loading balance.
		if (count % numprocs == rank) {
			int nodeId;
			istringstream lineStream(lineStr);
			while (lineStream >> nodeId) {
				rankMapCommunities.insert(make_pair(nodeId, count));
			}
		}

		++count;
	}

	communityFp.close();
}

template<class T>
void Reader<T>::getRankMapCommunity(unordered_map<T, int>& rankMapCommunities,
		unordered_map<int, int>& communitySizes, int rank, int numprocs) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getRankMapCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return;
	}

	// Community id starts from 0
	int count = 0;
	// Read until all the lines are read
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		// Each rank read its own communities. This method could also achieve rough loading balance.
		if (count % numprocs == rank) {
			int nodeId;
			int comSize = 0;
			istringstream lineStream(lineStr);
			while (lineStream >> nodeId) {
				rankMapCommunities.insert(make_pair(nodeId, count));
				++comSize;
			}

			communitySizes.insert(make_pair(count, comSize));
		}

		++count;
	}

	communityFp.close();
}

template<class T>
void Reader<T>::getDetectedMapCommunity(
		const unordered_map<T, int>& realMapCommunities,
		unordered_map<T, int>& disMapCommunities) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getDetectedMapCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return;
	}

	// Community id starts from 0
	int count = 0;
	// Read until all the lines are read
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		// Only read nodes that in the real communities
		int nodeId;
		typename unordered_map<T, T>::const_iterator it;
		istringstream lineStream(lineStr);
		while (lineStream >> nodeId) {
			it = realMapCommunities.find(nodeId);
			if (it != realMapCommunities.end()) {
				disMapCommunities.insert(make_pair(nodeId, count));
			}
		}

		++count;
	}

	communityFp.close();
}

/**
 * Get the network that only contains nodes of the corresponding communities
 */
template<class T>
double Reader<T>::getCommunityNetwork(
		const unordered_map<T, int>& disMapCommunities,
		unordered_map<T, unordered_map<T, double> >& communityNetwork,
		unordered_set<T>& outCommunityNodes, bool isUnweighted,
		bool isUndirected) {
	assert(this->fileName != "");
	ifstream networkFp;
	networkFp.open(this->fileName.c_str(), ios_base::in);
	if (!networkFp.is_open()) {
		fprintf(stderr, "getCommunityNetwork->cannot open %s\n", this->fileName.c_str());
		return -1;
	}

	double comNetWeight = 0;
	string lineStr;
	T srcId;
	T dstId;
	double weight;
	typedef unordered_map<T, double> innerMap;
	typename unordered_map<T, int>::const_iterator srcIt;
	typename unordered_map<T, int>::const_iterator dstIt;
	typename unordered_map<T, double>::iterator nbIter;
	while (networkFp.good()) {
		getline(networkFp, lineStr);
		if (lineStr.empty()) {
			continue;
		}

		istringstream lineStream(lineStr);
		lineStream >> srcId;
		lineStream >> dstId;

		// Ignore self-loop edge
		if (srcId == dstId) {
			continue;
		}
		weight = 1;

		// If weighted network, read in weight
		if (!isUnweighted) {
			lineStream >> weight;
		}

		// cout << srcId << "\t" << dstId << "\t" << weight << endl;

		srcIt = disMapCommunities.find(srcId);
		dstIt = disMapCommunities.find(dstId);
		if (srcIt != disMapCommunities.end()) {
			innerMap &nbs = communityNetwork[srcId];
			nbIter = nbs.find(dstId);

			if (nbIter == nbs.end()) {
				nbs.insert(make_pair(dstId, weight));
				comNetWeight += weight;
			}

			if (dstIt == disMapCommunities.end()) {
				outCommunityNodes.insert(dstId);
			} else {
				innerMap dstNbs = communityNetwork[dstId];
			}
		}

		// For undirected network, add the reverse edge
		if (isUndirected) {
			if (dstIt != disMapCommunities.end()) {
				innerMap &nbs = communityNetwork[dstId];
				nbIter = nbs.find(srcId);

				if (nbIter == nbs.end()) {
					nbs.insert(make_pair(srcId, weight));
					comNetWeight += weight;
				}

				if (srcIt == disMapCommunities.end()) {
					outCommunityNodes.insert(srcId);
				}
			}
		}
	}

	return comNetWeight;
}

/**
 * Get the network that only contains nodes of the corresponding communities
 */
template<class T>
double Reader<T>::getCommunityReversedNetwork(
		const unordered_map<T, int>& disMapCommunities,
		unordered_map<T, unordered_map<T, double> >& communityInNetwork,
		bool isUnweighted, bool isUndirected) {
	assert(this->fileName != "");
	ifstream networkFp;
	networkFp.open(this->fileName.c_str(), ios_base::in);
	if (!networkFp.is_open()) {
		fprintf(stderr, "getCommunityNetwork->cannot open %s\n", this->fileName.c_str());
		return -1;
	}

	double comNetWeight = 0;
	string lineStr;
	T srcId;
	T dstId;
	double weight;
	typedef unordered_map<T, double> innerMap;
	typename unordered_map<T, int>::const_iterator srcIt;
	typename unordered_map<T, int>::const_iterator dstIt;
	typename unordered_map<T, double>::iterator nbIter;
	while (networkFp.good()) {
		getline(networkFp, lineStr);
		if (lineStr.empty()) {
			continue;
		}

		istringstream lineStream(lineStr);
		lineStream >> srcId;
		lineStream >> dstId;

		// Ignore self-loop edge
		if (srcId == dstId) {
			continue;
		}
		weight = 1;

		// If weighted network, read in weight
		if (!isUnweighted) {
			lineStream >> weight;
		}

		// cout << srcId << "\t" << dstId << "\t" << weight << endl;

		srcIt = disMapCommunities.find(srcId);
		dstIt = disMapCommunities.find(dstId);
		if (dstIt != disMapCommunities.end()) {
			innerMap &nbs = communityInNetwork[dstId];
			nbIter = nbs.find(srcId);

			if (nbIter == nbs.end()) {
				nbs.insert(make_pair(srcId, weight));
				comNetWeight += weight;
			}

			if (srcIt != disMapCommunities.end()) {
				innerMap srcNbs = communityInNetwork[srcId];
			}
		}

		// For undirected network, add the reverse edge
		if (isUndirected) {
			if (srcIt != disMapCommunities.end()) {
				innerMap &nbs = communityInNetwork[srcId];
				nbIter = nbs.find(dstId);

				if (nbIter == nbs.end()) {
					nbs.insert(make_pair(dstId, weight));
					comNetWeight += weight;
				}
			}
		}
	}

	return comNetWeight;
}

template<class T>
void Reader<T>::getOutCommunityNodeInfo(unordered_map<int, int>& communitySizes,
		unordered_map<T, int>& disMapCommunities,
		const unordered_set<T>& outCommunityNodes) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getRankMapCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return;
	}

	// Community id starts from 0
	int count = 0;
	typename unordered_set<T>::const_iterator it;
	// Read until all the lines are read
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		int nodeId;
		bool flag = false;
		int comSize = 0;
		istringstream lineStream(lineStr);
		while (lineStream >> nodeId) {
			it = outCommunityNodes.find(nodeId);
			if (it != outCommunityNodes.end()) {
				flag = true;
				disMapCommunities.insert(make_pair(nodeId, count));
			}
			++comSize;
		}

		if (flag) {
			communitySizes.insert(make_pair(count, comSize));
		}

		++count;
	}

	communityFp.close();
}

/**
 * Get a community from a line of string
 */
template<class T>
unordered_set<T> Reader<T>::getLineCommunity(string lineStr) {
	unordered_set<T> community;
	int nodeId;
	istringstream lineStream(lineStr);
	while (lineStream >> nodeId) {
		community.insert(nodeId);
	}

	return community;
}

/**
 * Get the whole network
 */
template<class T>
double Reader<T>::getNetwork(bool isUnweighted, bool isUndirected,
		unordered_map<T, unordered_map<T, double> >& network) {
	assert(this->fileName != "");
	ifstream networkFp;
	networkFp.open(this->fileName.c_str(), ios_base::in);
	if (!networkFp.is_open()) {
		fprintf(stderr, "getNetwork->cannot open %s\n", this->fileName.c_str());
		return -1;
	}

	double totalWeight = 0;
	string lineStr;
	T srcId;
	T dstId;
	double weight;
	typedef unordered_map<T, double> innerMap;
	typename unordered_map<T, double>::iterator nbIter;
	while (networkFp.good()) {
		getline(networkFp, lineStr);
		if (lineStr.empty()) {
			continue;
		}

		istringstream lineStream(lineStr);
		lineStream >> srcId;
		lineStream >> dstId;

		// Ignore self-loop edge
		if (srcId == dstId) {
			continue;
		}
		weight = 1;

		// If weighted network, read in weight
		if (!isUnweighted) {
			lineStream >> weight;
		}

		// cout << srcId << "\t" << dstId << "\t" << weight << endl;

		innerMap &nbs = network[srcId];
		nbIter = nbs.find(dstId);

		if (nbIter == nbs.end()) {
			nbs.insert(make_pair(dstId, weight));
			totalWeight += weight;
		}

		innerMap dstNbs = network[dstId];

		// For undirected network, add the reverse edge
		if (isUndirected) {
			innerMap &nbs = network[dstId];
			nbIter = nbs.find(srcId);

			if (nbIter == nbs.end()) {
				nbs.insert(make_pair(srcId, weight));
				totalWeight += weight;
			}
		}
	}

	return totalWeight;
}

/**
 * Get the whole network
 */
template<class T>
double Reader<T>::getReversedNetwork(bool isUnweighted, bool isUndirected,
		unordered_map<T, unordered_map<T, double> >& inNet) {
	assert(this->fileName != "");
	ifstream networkFp;
	networkFp.open(this->fileName.c_str(), ios_base::in);
	if (!networkFp.is_open()) {
		fprintf(stderr, "getNetwork->cannot open %s\n", this->fileName.c_str());
		return -1;
	}

	double totalWeight = 0;
	string lineStr;
	T srcId;
	T dstId;
	double weight;
	typedef unordered_map<T, double> innerMap;
	typename unordered_map<T, double>::iterator nbIter;
	while (networkFp.good()) {
		getline(networkFp, lineStr);
		if (lineStr.empty()) {
			continue;
		}

		istringstream lineStream(lineStr);
		lineStream >> srcId;
		lineStream >> dstId;

		// Ignore self-loop edge
		if (srcId == dstId) {
			continue;
		}
		weight = 1;

		// If weighted network, read in weight
		if (!isUnweighted) {
			lineStream >> weight;
		}

		// cout << srcId << "\t" << dstId << "\t" << weight << endl;

		innerMap &nbs = inNet[dstId];
		nbIter = nbs.find(srcId);

		if (nbIter == nbs.end()) {
			nbs.insert(make_pair(srcId, weight));
			totalWeight += weight;
		}

		innerMap srcNbs = inNet[srcId];

		// For undirected network, add the reverse edge
		if (isUndirected) {
			innerMap &nbs = inNet[srcId];
			nbIter = nbs.find(dstId);

			if (nbIter == nbs.end()) {
				nbs.insert(make_pair(dstId, weight));
				totalWeight += weight;
			}
		}
	}

	return totalWeight;
}

/**
 * Get the whole community and stored in vector
 */
template<class T>
long Reader<T>::getCommunity(vector<unordered_set<T> >& communities) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return -1;
	}

	// Read until all the lines are read
	long numNodes = 0;
	int count = 0;
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		unordered_set<T> community = this->getLineCommunity(lineStr);
		communities.push_back(community);
		numNodes += community.size();
		++count;
	}

	communityFp.close();

	return numNodes;
}

/**
 * Get the community of nodes
 */
template<class T>
long Reader<T>::getMapCommunity(unordered_map<T, int>& communities) {
	assert(this->fileName != "");
	ifstream communityFp;
	string lineStr;
	communityFp.open(this->fileName.c_str(), ios_base::in);
	if (!communityFp.is_open()) {
		fprintf(stderr, "getMapCommunity->cannot open %s\n", this->fileName.c_str());
		//cerr << "getRankCommunity->cannot open the file!" << this->fileName;
		return -1;
	}

	// Community id starts from 0
	long numNodes = 0;
	int count = 0;
	// Read until all the lines are read
	while (communityFp.good()) {
		getline(communityFp, lineStr);
		if (lineStr.empty()) {
			//cout << "linStr is empty" << endl;
			continue;
		}

		int nodeId;
		istringstream lineStream(lineStr);
		while (lineStream >> nodeId) {
			communities.insert(make_pair(nodeId, count));
			++numNodes;
		}

		++count;
	}

	communityFp.close();
	return numNodes;
}

template<class T>
Reader<T>::~Reader(void) {
}

#endif
