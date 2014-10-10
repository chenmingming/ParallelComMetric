/*
 * MPIMetric.cpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Miles
 */

#ifndef MPIMETRIC_CPP_
#define MPIMETRIC_CPP_

#include "MPIMetric.h"

template<class T>
MPIMetric<T>::MPIMetric() {

}

/**
 * Calculate metrics without ground truth communities. The metrics include modularity,
 * Modularity Density, intra-edges, intra-density, contraction, inter-edges, expansion, and conductance
 */
template<class T>
void MPIMetric<T>::computeMetricWithoutGroundTruth(int rankid,
		double comNetWeight, unordered_map<int, int>& communitySizes,
		unordered_map<T, int>& disMapCommunities,
		const unordered_map<T, unordered_map<T, double> >& communityNetwork,
		const unordered_map<T, unordered_map<T, double> >& communityInNetwork,
		bool isUndirected, double *metrics) {
	double totalWeight = 0;
	MPI_Allreduce(&comNetWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	unordered_map<int, unordered_map<int, double> > communityWeights;
	unordered_map<int, unordered_map<int, double> > communityEdges;
	typename unordered_map<T, unordered_map<T, double> >::const_iterator netIt;

	for (netIt = communityNetwork.begin(); netIt != communityNetwork.end();
			++netIt) {
		int srcId = netIt->first;
		int srcComId = disMapCommunities[srcId];
		// The weights between a certain community to its neighboring communities
		unordered_map<int, double> &comWeights = communityWeights[srcComId];
		unordered_map<int, double> &comEdges = communityEdges[srcComId];
		unordered_map<T, double> nbs = netIt->second;
		typename unordered_map<T, double>::const_iterator nbIt;
		for (nbIt = nbs.begin(); nbIt != nbs.end(); ++nbIt) {
			int dstId = nbIt->first;
			int dstComId = disMapCommunities[dstId];
			double weight = nbIt->second;
			double comWeight = comWeights[dstComId];
			comWeights[dstComId] = comWeight + weight;
			double comEdge = comEdges[dstComId];
			comEdges[dstComId] = comEdges[dstComId] + 1;

#ifdef MPI_DEBUG
			if (OUTPUT_RANK == rankid) {
				cout << "comWeight = " << comWeight
				<< ", comWeights[dstComId] = " << comWeights[dstComId]
				<< endl;
			}
#endif
		}
	}

	// directed network
	unordered_map<int, double> communityInWeights;
	if (!isUndirected) {
		for (netIt = communityInNetwork.begin();
				netIt != communityInNetwork.end(); ++netIt) {
			int srcId = netIt->first;
			int srcComId = disMapCommunities[srcId];
			unordered_map<T, double> nbs = netIt->second;
			typename unordered_map<T, double>::const_iterator nbIt;
			for (nbIt = nbs.begin(); nbIt != nbs.end(); ++nbIt) {
				int dstId = nbIt->first;
				unsigned int isExist = disMapCommunities.count(dstId);
				bool outComNodes = false;
				if (isExist) {
					int dstComId = disMapCommunities[dstId];
					if (srcComId != dstComId) {
						outComNodes = true;
					}
				} else {
					outComNodes = true;
				}

				if (outComNodes) {
					double weight = nbIt->second;
					double comWeight = communityInWeights[srcComId];
					communityInWeights[srcComId] = communityInWeights[srcComId]
							+ weight;
				}

			}
		}
	}

	double rankModularity = 0;
	double rankQds = 0;
	double rankIntraEdges = 0;
	double rankIntraDensity = 0;
	double rankContraction = 0;
	double rankInterEdges = 0;
	//double interDensity = 0;
	double rankExpansion = 0;
	double rankConductance = 0;

	// Calculate metrics for this rank
	typename unordered_map<int, unordered_map<int, double> >::const_iterator comIter;
	for (comIter = communityWeights.begin(); comIter != communityWeights.end();
			++comIter) {
		int srcComId = comIter->first;
		unordered_map<int, double> comWeights = comIter->second;
		unordered_map<int, double> comEdges = communityEdges[srcComId];

#ifdef MPI_DEBUG
		if (OUTPUT_RANK == rankid) {
			cout << "rankComSize = " << communityWeights.size()
			<< ", srcComId = " << srcComId << ", nbComNum = "
			<< comWeights.size() << endl;
		}
#endif

		double inWeights = comWeights[srcComId];
		double inEdges = comEdges[srcComId];
		double srcComSize = communitySizes[srcComId];
		double outWeights = 0;
		double out_incoming_weights = communityInWeights[srcComId];
		double splitPenalty = 0;
		typename unordered_map<int, double>::const_iterator weightIter;
		for (weightIter = comWeights.begin(); weightIter != comWeights.end();
				++weightIter) {
			int dstComId = weightIter->first;

			if (srcComId != dstComId) {
				int dstComSize = communitySizes[dstComId];
				double comWeight = weightIter->second;
				outWeights += comWeight;
				double comEdge = comEdges[dstComId];
				double sp = (comWeight / totalWeight)
						* (comEdge / (srcComSize * dstComSize));
				splitPenalty += sp;
			}
		}

//		if (OUTPUT_RANK == rankid) {
//			cout << "inWeights = " << inWeights << ", outWeights = "
//					<< outWeights << ", out_incoming_weights = "
//					<< out_incoming_weights << endl;
//		}

// modularity
		if (isUndirected) {
			rankModularity += inWeights / totalWeight
					- pow((inWeights + outWeights) / totalWeight, 2);
		} else {
			rankModularity += inWeights / totalWeight
					- ((inWeights + outWeights)
							* (inWeights + out_incoming_weights))
							/ pow(totalWeight, 2);
		}

		// Modularity Density Qds
		double inDensity = 0;
		if (srcComSize > 1) {
			inDensity = inEdges / (srcComSize * (srcComSize - 1));
		}
		if (isUndirected) {
			rankQds += (inWeights / totalWeight) * inDensity
					- pow(((inWeights + outWeights) / totalWeight) * inDensity,
							2) - splitPenalty;
		} else {
			rankQds += (inWeights / totalWeight) * inDensity
					- (((inWeights + outWeights)
							* (inWeights + out_incoming_weights))
							/ pow(totalWeight, 2)) * pow(inDensity, 2)
					- splitPenalty;
		}

		// intra-edges
		if (isUndirected) {
			rankIntraEdges += inWeights / 2;
		} else {
			rankIntraEdges += inWeights;
		}
		// contraction: average degree
		if (inWeights == 0 || srcComSize == 0) {
			rankContraction += 0;
		} else {
			rankContraction += inWeights / srcComSize;
		}
		// intra-density
		rankIntraDensity += inDensity;

		rankInterEdges += outWeights;
		// inter-density
		// if (numNodes == srcComSize) {
		// interDensity += 0;
		// } else {
		// interDensity += outWeights
		// / (srcComSize * (numNodes - srcComSize));
		// }
		if (outWeights == 0 || srcComSize == 0) {
			rankExpansion += 0;
		} else {
			rankExpansion += outWeights / srcComSize;
		}

		// Avoid that totalInterEdges==0 and communityEdges[i][i]==0
		if (outWeights == 0) {
			rankConductance += 0;
		} else {
			rankConductance += outWeights / (inWeights + outWeights);
		}
	} // community for

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankQ = " << rankModularity << ", rankQds = " << rankQds
		<< ", rankIntraEdges = " << rankIntraEdges
		<< ", rankIntraDensity = " << rankIntraDensity
		<< ", rankContraction = " << rankContraction
		<< ", rankInterEdges = " << rankInterEdges
		<< ", rankExpansion = " << rankExpansion
		<< ", rankConductance = " << rankConductance << endl;
	}
#endif

	// Get the total values of the metrics
	int rankComNum = communityWeights.size();
	int numComs = 0;
	double modularity = 0;
	double Qds;
	double intraEdges = 0;
	double intraDensity = 0;
	double contraction = 0;
	double interEdges = 0;
	double expansion = 0;
	double conductance = 0;

	MPI_Allreduce(&rankComNum, &numComs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankModularity, &modularity, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankQds, &Qds, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankIntraEdges, &intraEdges, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankIntraDensity, &intraDensity, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankContraction, &contraction, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankInterEdges, &interEdges, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankExpansion, &expansion, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankConductance, &conductance, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	metrics[0] = modularity;
	metrics[1] = Qds;
	metrics[2] = intraEdges / numComs;
	metrics[3] = intraDensity / numComs;
	metrics[4] = contraction / numComs;
	metrics[5] = interEdges / numComs;
	metrics[6] = expansion / numComs;
	metrics[7] = conductance / numComs;

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "Q = " << metrics[0] << ", Qds = " << metrics[1]
		<< ", intraEdges = " << metrics[2] << ", intraDensity = "
		<< metrics[3] << ", contraction = " << metrics[4]
		<< ", interEdges = " << metrics[5] << ", expansion = "
		<< metrics[6] << ", conductance = " << metrics[7] << endl;
	}
#endif
}

/**
 * Compute Variation Information (VI) and Normalized Mutual Information
 * (NMI): Disjoint community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realCommunities: the real communities of this rank
 * @param disCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (VI and NMI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeInfoEntropyMetric(int rankid, int numprocs,
		const vector<unordered_set<T> >& realCommunities,
		const vector<unordered_set<T> >& disCommunities, double *metrics,
		double *entropyCycles) {
	int realComNum = realCommunities.size();
	int disComNum = disCommunities.size();

// Get the number of nodes in this rank
	double rankNumNodes = 0;
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> community = realCommunities[i];
		rankNumNodes += community.size();
	}

// Get the total number of nodes of the network
	double numNodes = 0;
	MPI_Allreduce(&rankNumNodes, &numNodes, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "Into computeInfoEntropyMetric with realComNum = " << realComNum
		<< " and disComNum = " << disComNum << ", rankNumNodes = "
		<< rankNumNodes << ", numNodes = " << numNodes << endl;
	}
#endif

// Compute x entropy or entropy of ground truth communities
	double startCompCycle = rdtsc();
	double rankXentropy = 0;
	for (int i = 0; i < realComNum; ++i) {
		double ni = realCommunities[i].size();
		rankXentropy += -(ni / numNodes) * log2(ni / numNodes);
	}

// Compute y entropy or entropy of detected communities and also get the
// number of items should be sent of this rank
	unsigned int num_item_sent = disComNum;
	double rankYentropy = 0;
	for (int j = 0; j < disComNum; ++j) {
		double nj = disCommunities[j].size();
		num_item_sent += nj;
		rankYentropy += -(nj / numNodes) * log2(nj / numNodes);
	}

// Compute VI and NMI for this rank.
	double rankVI = 0;
	double rankNMI = 0;
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> realCommunity = realCommunities[i];
		double ni = realCommunity.size();

		for (int j = 0; j < disComNum; ++j) {
			unordered_set<T> disCommunity = disCommunities[j];
			double nj = disCommunity.size();
			double nij = this->getIntersectionNum(realCommunity, disCommunity);

			if (nij != 0) {
				rankVI += nij * log2((nij * nij) / (ni * nj));
				rankNMI += (nij / numNodes)
						* log2((nij * numNodes) / (ni * nj));
			}
		}
	}
	double endCompCycle = rdtsc();
	entropyCycles[0] += endCompCycle - startCompCycle;

// If numprocs<=1, do not conduct the message sending and receiving part
	if (numprocs <= 1) {
		double VI = -rankVI / numNodes;
		double NMI = (2 * rankNMI) / (rankXentropy + rankYentropy);
		metrics[0] = VI;
		metrics[1] = NMI;
		return;
	}

	double startMsgCycle, endMsgCycle;
	unsigned int message_send_num = 0;
	int send_flag, recv_flag = 0;
	MPI_Status status;

	/* Message send variables */
	T *send_buff = (T *) malloc(num_item_sent * sizeof(T));

// Initialize the send_buff to be the detected communities of this rank
	int step = 0;
	for (int i = 0; i < disComNum; ++i) {
		unordered_set<T> community = disCommunities[i];
		typename unordered_set<T>::const_iterator comIter;
		for (comIter = community.begin(); comIter != community.end();
				++comIter) {
			send_buff[step++] = *comIter;
		}

		// Use -1 as the boundary of two communities
		send_buff[step++] = -1;
	}

#ifdef MPI_DEBUG
// Output self send_buff
	if (rankid != OUTPUT_RANK) {
		cout << "self send_buff: " << endl;
		for (int i = 0; i < num_item_sent; ++i) {
			if (send_buff[i] == -1) {
				cout << endl;
			} else {
				cout << send_buff[i] << " ";
			}
		}
	}
#endif

	MPI_Request send_request = MPI_REQUEST_NULL;
	int send_dest = (rankid + numprocs - 1) % numprocs;

	/* Message recv variables */
	T *recv_buff;
	int num_item_recvd = 0;
	MPI_Request recv_request = MPI_REQUEST_NULL;
	int recv_source = (rankid + 1) % numprocs;

	/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
	 * next MPI rank at the first step and also the message tag is its rankid.
	 */
	int recv_tag = rankid;

// Save the received communities
	vector<unordered_set<T> > bufCommunities;

	/* When receiving (numprocs - 1) messages, exit */
	while (message_send_num < (numprocs - 1)) {
		startMsgCycle = rdtsc();
		/* At first step, send its own slice of matrix B */
		/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
		if (recv_tag != send_dest) {
			MPI_Isend(send_buff, num_item_sent, MPI_INT, send_dest, recv_tag,
					MPI_COMM_WORLD, &send_request);
		}

		recv_flag = 0;
		while (!recv_flag) {
			MPI_Iprobe(recv_source, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_flag,
					&status);

			// When ready to receive
			if (recv_flag) {
				// Get the number of item that should be received
				MPI_Get_count(&status, MPI_INT, &num_item_recvd);
				// Post a recv to receive the message
				recv_buff = (T *) malloc(num_item_recvd * sizeof(T));
				MPI_Recv(recv_buff, num_item_recvd, MPI_INT, recv_source,
						MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				/* recv_tag is the taskid that owns this slice of matrix B */
				recv_tag = status.MPI_TAG;
			}
		}

		/* wait until the MPI_IRecv and MPI_ISend succeed */
		send_flag = 0;
		while (!send_flag) {
			/* Only when send_request is not null, shall we test its status */
			if (send_request != MPI_REQUEST_NULL) {
				MPI_Test(&send_request, &send_flag, &status);

				/* When finishing sending message, post a recv request */
				if (send_flag) {
					++message_send_num;
					free(send_buff);
				}
			} /* if request */
		} /* while flag */
		endMsgCycle = rdtsc();
		entropyCycles[1] += endMsgCycle - startMsgCycle;

#ifdef MPI_DEBUG
		// Output recv_buff
		if (OUTPUT_RANK == rankid) {
			cout << "recv_buff: " << endl;
			for (int i = 0; i < num_item_recvd; ++i) {
				if (recv_buff[i] == -1) {
					cout << endl;
				} else {
					cout << recv_buff[i] << " ";
				}
			}
		}
#endif

		/* Convert recv_buff to matrix vector<unordered_set<T> > disCommunities and copy to send_buff */
		if (message_send_num < (numprocs - 1) && recv_tag != send_dest) {
			num_item_sent = num_item_recvd;
			send_buff = (T *) malloc(num_item_sent * sizeof(T));
		}
		bufCommunities.clear();
		unordered_set<T> community;
		for (int i = 0; i < num_item_recvd; ++i) {
			if (recv_buff[i] == -1) {
				bufCommunities.push_back(community);
				community.clear();
			} else {
				community.insert(recv_buff[i]);
			}

			if (message_send_num < (numprocs - 1) && recv_tag != send_dest) {
				send_buff[i] = recv_buff[i];
			}
		}

		free(recv_buff);

#ifdef MPI_DEBUG
		// Output received communities
		if (OUTPUT_RANK == rankid) {
			int disComSize = bufCommunities.size();
			cout << disComSize << " discovered communities: " << endl;
			for (int i = 0; i < disComSize; ++i) {
				unordered_set<int> community = bufCommunities[i];
				typename unordered_set<int>::const_iterator comIter;
				for (comIter = community.begin(); comIter != community.end();
						++comIter) {
					cout << *comIter << " ";
				}
				cout << endl;
			}
		}
#endif

		// Compute the VI and NMI from other ranks
		startCompCycle = rdtsc();
		disComNum = bufCommunities.size();
		for (int i = 0; i < realComNum; ++i) {
			unordered_set<T> realCommunity = realCommunities[i];
			double ni = realCommunity.size();

			for (int j = 0; j < disComNum; ++j) {
				unordered_set<T> disCommunity = bufCommunities[j];
				double nj = disCommunity.size();
				double nij = this->getIntersectionNum(realCommunity,
						disCommunity);

				if (nij != 0) {
					rankVI += nij * log2((nij * nij) / (ni * nj));
					rankNMI += (nij / numNodes)
							* log2((nij * numNodes) / (ni * nj));
				}
			}
		}
		endCompCycle = rdtsc();
		entropyCycles[0] += endCompCycle - startCompCycle;
	} /* while */

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankVI = " << rankVI << ", rankNMI = " << rankNMI << endl;
	}
#endif

	double VI = 0;
	double xentropy = 0;
	double yentropy = 0;
	double NMI = 0;
	MPI_Allreduce(&rankVI, &VI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankXentropy, &xentropy, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankYentropy, &yentropy, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankNMI, &NMI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	VI = -VI / numNodes;
	NMI = (2 * NMI) / (xentropy + yentropy);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "VI = " << VI << ", xentropy = " << xentropy << ", yentropy = "
		<< yentropy << ", NMI = " << NMI << endl;
	}
#endif

	metrics[0] = VI;
	metrics[1] = NMI;
}

/**
 * Compute Variation Information (VI) and Normalized Mutual Information
 * (NMI): Disjoint community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realCommunities: the real communities of this rank
 * @param disCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (VI and NMI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeInfoEntropyMetric_backup(int rankid, int numprocs,
		const vector<unordered_set<T> >& realCommunities,
		const vector<unordered_set<T> >& disCommunities, double *metrics) {
	int realComNum = realCommunities.size();
	int disComNum = disCommunities.size();

// Get the number of nodes in this rank
	double rankNumNodes = 0;
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> community = realCommunities[i];
		rankNumNodes += community.size();
	}

// Get the total number of nodes of the network
	double numNodes = 0;
	MPI_Allreduce(&rankNumNodes, &numNodes, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "Into computeInfoEntropyMetric with realComNum = " << realComNum
		<< " and disComNum = " << disComNum << ", rankNumNodes = "
		<< rankNumNodes << ", numNodes = " << numNodes << endl;
	}
#endif

// Compute x entropy or entropy of ground truth communities
	double rankXentropy = 0;
	for (int i = 0; i < realComNum; ++i) {
		double ni = realCommunities[i].size();
		rankXentropy += -(ni / numNodes) * log2(ni / numNodes);
	}

// Compute y entropy or entropy of detected communities and also get the
// max number of nodes of a community among all the communities of this rank
	int rankMaxComSize = -1;
	double rankYentropy = 0;
	for (int j = 0; j < disComNum; ++j) {
		double nj = disCommunities[j].size();
		if (nj > rankMaxComSize) {
			rankMaxComSize = nj;
		}
		rankYentropy += -(nj / numNodes) * log2(nj / numNodes);
	}

// Get the max number of nodes of a community among all the detected communities
	int maxComSize = 0;
	MPI_Allreduce(&rankMaxComSize, &maxComSize, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

// Get the max number of communities that a rank may have
	int maxComNum = 0;
	MPI_Allreduce(&disComNum, &maxComNum, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "maxComNum = " << maxComNum << ", maxComSize = " << maxComSize
		<< endl;
	}
#endif

// Compute VI and NMI for this rank.
	double rankVI = 0;
	double rankNMI = 0;
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> realCommunity = realCommunities[i];
		double ni = realCommunity.size();

		for (int j = 0; j < disComNum; ++j) {
			unordered_set<T> disCommunity = disCommunities[j];
			double nj = disCommunity.size();
			double nij = this->getIntersectionNum(realCommunity, disCommunity);

			if (nij != 0) {
				rankVI += nij * log2((nij * nij) / (ni * nj));
				rankNMI += (nij / numNodes)
						* log2((nij * numNodes) / (ni * nj));
			}
		}
	}

// If numprocs<=1, do not conduct the message sending and receiving part
	if (numprocs <= 1) {
		double VI = -rankVI / numNodes;
		double NMI = (2 * rankNMI) / (rankXentropy + rankYentropy);

		metrics[0] = VI;
		metrics[1] = NMI;
		return;
	}

	int message_send_num = 0;
	int send_flag, recv_flag = 0;
	MPI_Status status;

// The total number of elements will be sent or received
	unsigned int count = maxComNum * maxComSize;

	/* Message send variables */
//T **send_buff = malloc2d(maxComNum, maxComSize);
	T send_buff[maxComNum][maxComSize];

// Initialize the send_buff to be the detected communities of this rank
	for (int i = 0; i < disComNum; ++i) {
		unordered_set<T> community = disCommunities[i];
		int step = 0;
		typename unordered_set<T>::const_iterator comIter;
		for (comIter = community.begin(); comIter != community.end();
				++comIter) {
			send_buff[i][step] = *comIter;
			++step;
		}

		// Fill the remaining elements with END_OF_COMMUNITY
		for (int j = step; j < maxComSize; ++j) {
			send_buff[i][j] = END_OF_COMMUNITY;
		}
	}

// Fill the remaining rows with END_OF_COMMUNITY
	for (int i = disComNum; i < maxComNum; ++i) {
		for (int j = 0; j < maxComSize; ++j) {
			send_buff[i][j] = END_OF_COMMUNITY;
		}
	}

#ifdef MPI_DEBUG
// Output self send_buff
	if (rankid != OUTPUT_RANK) {
		cout << "self send_buff: " << endl;
		for (int i = 0; i < maxComNum; ++i) {
			for (int j = 0; j < maxComSize; ++j) {
				cout << send_buff[i][j] << " ";
			}
			cout << endl;
		}
	}
#endif

	MPI_Request send_request = MPI_REQUEST_NULL;
	int send_dest = (rankid + numprocs - 1) % numprocs;

	/* Message recv variables */
//T **recv_buff = malloc2d(maxComNum, maxComSize);
	T recv_buff[maxComNum][maxComSize];
	MPI_Request recv_request = MPI_REQUEST_NULL;
	int recv_source = (rankid + 1) % numprocs;

	/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
	 * next MPI rank at the first step and also the message tag is its rankid.
	 */
	int recv_tag = rankid;

// Save the received communities
	vector<unordered_set<T> > bufCommunities;

	/* When receiving (numprocs - 1) messages, exit */
	while (message_send_num < (numprocs - 1)) {
		/* First, post a recv request to prevent deadlock */
		MPI_Irecv(&recv_buff[0][0], count, MPI_INT, recv_source, MPI_ANY_TAG,
				MPI_COMM_WORLD, &recv_request);

		/* At first step, send its own slice of matrix B */
		/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
		if (recv_tag != send_dest) {
			MPI_Isend(&send_buff[0][0], count, MPI_INT, send_dest, recv_tag,
					MPI_COMM_WORLD, &send_request);
		}

		/* wait until the MPI_IRecv and MPI_ISend succeed */
		send_flag = 0;
		recv_flag = 0;
		while (!recv_flag || !send_flag) {
			/* Only when recv_request is not null, shall we test its status */
			if (recv_request != MPI_REQUEST_NULL) {
				MPI_Test(&recv_request, &recv_flag, &status);

				/* When finishing recving message, post a send request */
				if (recv_flag) {
					/* recv_tag is the taskid that owns this slice of matrix B */
					recv_tag = status.MPI_TAG;
				} /* if flag */
			} /* if request*/

			/* Only when send_request is not null, shall we test its status */
			if (send_request != MPI_REQUEST_NULL) {
				MPI_Test(&send_request, &send_flag, &status);

				/* When finishing sending message, post a recv request */
				if (send_flag) {
					++message_send_num;
				}
			} /* if request */
		} /* while flag */

#ifdef MPI_DEBUG
		// Output recv_buff
		if (OUTPUT_RANK == rankid) {
			cout << "recv_buff: " << endl;
			for (int i = 0; i < maxComNum; ++i) {
				for (int j = 0; j < maxComSize; ++j) {
					cout << recv_buff[i][j] << " ";
				}
				cout << endl;
			}
		}
#endif

		/* Convert recv_buff to matrix vector<unordered_set<T> > disCommunities */
		bufCommunities.clear();
		for (int i = 0; i < maxComNum; ++i) {
			// If this line does not correspond to a community, ignore it
			if (recv_buff[i][0] == END_OF_COMMUNITY) {
				for (int j = 0; j < maxComSize; ++j) {
					send_buff[i][j] = recv_buff[i][j];
				}
				continue;
			}
			unordered_set<T> community;
			for (int j = 0; j < maxComSize; ++j) {
				int nodeId = recv_buff[i][j];
				send_buff[i][j] = nodeId;
				if (nodeId != END_OF_COMMUNITY) {
					community.insert(nodeId);
				}
			}

			bufCommunities.push_back(community);
		}

#ifdef MPI_DEBUG
		// Output received communities
		if (OUTPUT_RANK == rankid) {
			int disComSize = bufCommunities.size();
			cout << disComSize << " discovered communities: " << endl;
			for (int i = 0; i < disComSize; ++i) {
				unordered_set<int> community = bufCommunities[i];
				typename unordered_set<int>::const_iterator comIter;
				for (comIter = community.begin(); comIter != community.end();
						++comIter) {
					cout << *comIter << " ";
				}
				cout << endl;
			}
		}
#endif

		// Compute the VI and NMI from other ranks
		disComNum = bufCommunities.size();
		for (int i = 0; i < realComNum; ++i) {
			unordered_set<T> realCommunity = realCommunities[i];
			double ni = realCommunity.size();

			for (int j = 0; j < disComNum; ++j) {
				unordered_set<T> disCommunity = bufCommunities[j];
				double nj = disCommunity.size();
				double nij = this->getIntersectionNum(realCommunity,
						disCommunity);

				if (nij != 0) {
					rankVI += nij * log2((nij * nij) / (ni * nj));
					rankNMI += (nij / numNodes)
							* log2((nij * numNodes) / (ni * nj));
				}
			}
		}
	} /* while */

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankVI = " << rankVI << ", rankNMI = " << rankNMI << endl;
	}
#endif

// free memory
//free2d((T **) send_buff);
//free2d((T **) recv_buff);

	double VI = 0;
	double xentropy = 0;
	double yentropy = 0;
	double NMI = 0;
	MPI_Allreduce(&rankVI, &VI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankXentropy, &xentropy, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankYentropy, &yentropy, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	MPI_Allreduce(&rankNMI, &NMI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	VI = -VI / numNodes;
	NMI = (2 * NMI) / (xentropy + yentropy);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "VI = " << VI << ", xentropy = " << xentropy << ", yentropy = "
		<< yentropy << ", NMI = " << NMI << endl;
	}
#endif

	metrics[0] = VI;
	metrics[1] = NMI;
}

//----------------------------------------------------------------------------------------//
/**
 * Compute F-measure and Normalized Van Dongen metric (NVD): Disjoint
 * community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realCommunities: the real communities of this rank
 * @param disCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (VI and NMI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeClusterMatchingMetric(int rankid, int numprocs,
		const vector<unordered_set<T> >& realCommunities,
		const vector<unordered_set<T> >& disCommunities, double *metrics,
		double *clusteringCycles) {
	int realComNum = realCommunities.size();
	int disComNum = disCommunities.size();

// Get the number of nodes in this rank and also the max number of nodes among real communities
	double rankNumNodes = 0;
	int real_num_item_sent = realComNum;
	for (int i = 0; i < realComNum; ++i) {
		int ni = realCommunities[i].size();
		rankNumNodes += ni;
		real_num_item_sent += ni;
	}

// Get the max number of nodes of a community among all the detected communities of this rank
	int dis_num_item_sent = disComNum;
	for (int j = 0; j < disComNum; ++j) {
		int nj = disCommunities[j].size();
		dis_num_item_sent += nj;
	}

// Get the total number of nodes of the network
	double numNodes = 0;
	MPI_Allreduce(&rankNumNodes, &numNodes, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

// Compute NVD and FMeasure for this rank.
	double rankNVD = 0;
	double rankFMeasure = 0;

// (1) To get FMeasure and the first part of NVD
// To record the max of each real community of this rank
	double startCompCycle = rdtsc();

	//double realMaxFMeasures[realComNum];
	double *realMaxFMeasures = (double *) malloc(realComNum * sizeof(double));
//	double realMaxCommons[realComNum];
	double *realMaxCommons = (double *) malloc(realComNum * sizeof(double));
	memset(realMaxFMeasures, 0, realComNum * sizeof(double));
	memset(realMaxCommons, 0, realComNum * sizeof(double));
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> realCommunity = realCommunities[i];
		double ni = realCommunity.size();
		double maxCommon = -1;
		double maxFMeasure = -1;

		for (int j = 0; j < disComNum; ++j) {
			unordered_set<T> disCommunity = disCommunities[j];
			double nj = disCommunity.size();
			double nij = this->getIntersectionNum(realCommunity, disCommunity);

			if (nij > maxCommon) {
				maxCommon = nij;
			}

			double tmpFMeasure = (2 * nij) / (ni + nj);
			if (tmpFMeasure > maxFMeasure) {
				maxFMeasure = tmpFMeasure;
			}
		}

		if (maxFMeasure > realMaxFMeasures[i]) {
			realMaxFMeasures[i] = maxFMeasure;
		}

		if (maxCommon > realMaxCommons[i]) {
			realMaxCommons[i] = maxCommon;
		}
//		rankNVD += maxCommon;
//		rankFMeasure += ni * maxFMeasure;
	}
	double endCompCycle = rdtsc();
	clusteringCycles[0] += endCompCycle - startCompCycle;

	int message_send_num = 0;
	double startMsgCycle, endMsgCycle;
	int send_flag, recv_flag = 0;
	MPI_Status status;

	/* Message send and recv variables */
	int dis_num_item_recvd = 0;
	int recv_source = (rankid + 1) % numprocs;

	/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
	 * next MPI rank at the first step and also the message tag is its rankid.
	 */
	int recv_tag = rankid;

// Save the received communities
	vector<unordered_set<T> > bufCommunities;

	/* Message send variables */
	MPI_Request send_request = MPI_REQUEST_NULL;
	int send_dest = (rankid + numprocs - 1) % numprocs;

// Only if numprocs>1, shall we consider the message passing part
	if (numprocs > 1) {
		T *recv_buff;
		T *send_buff = (T *) malloc(dis_num_item_sent * sizeof(T));

		// Initialize the send_buff to be the detected communities of this rank
		int step = 0;
		for (int i = 0; i < disComNum; ++i) {
			unordered_set<T> community = disCommunities[i];
			typename unordered_set<T>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				send_buff[step++] = *comIter;
			}

			send_buff[step++] = -1;
		}

#ifdef MPI_DEBUG
		// Output self send_buff
		if (rankid == OUTPUT_RANK) {
			cout << "self send_buff: " << endl;
			for (int i = 0; i < dis_num_item_sent; ++i) {
				if (send_buff[i] == -1) {
					cout << endl;
				} else {
					cout << send_buff[i] << " ";
				}
			}
		}
#endif

		/* When receiving (numprocs - 1) messages, exit */
		while (message_send_num < (numprocs - 1)) {
			startMsgCycle = rdtsc();
			/* At first step, send its own slice of matrix B */
			/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
			if (recv_tag != send_dest) {
				MPI_Isend(send_buff, dis_num_item_sent, MPI_INT, send_dest,
						recv_tag, MPI_COMM_WORLD, &send_request);
			}

			recv_flag = 0;
			while (!recv_flag) {
				MPI_Iprobe(recv_source, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_flag,
						&status);

				// When ready to receive
				if (recv_flag) {
					// Get the number of item that should be received
					MPI_Get_count(&status, MPI_INT, &dis_num_item_recvd);
					// Post a recv to receive the message
					recv_buff = (T *) malloc(dis_num_item_recvd * sizeof(T));
					MPI_Recv(recv_buff, dis_num_item_recvd, MPI_INT,
							recv_source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					/* recv_tag is the taskid that owns this slice of matrix B */
					recv_tag = status.MPI_TAG;
				}
			}

			/* wait until the MPI_IRecv and MPI_ISend succeed */
			send_flag = 0;
			while (!send_flag) {
				/* Only when send_request is not null, shall we test its status */
				if (send_request != MPI_REQUEST_NULL) {
					MPI_Test(&send_request, &send_flag, &status);

					/* When finishing sending message, post a recv request */
					if (send_flag) {
						++message_send_num;
						free(send_buff);

//						cout << "First: message_send_num = " << message_send_num
//								<< endl;
					}
				} /* if request */
			} /* while flag */
			endMsgCycle = rdtsc();
			clusteringCycles[1] += endMsgCycle - startMsgCycle;

#ifdef MPI_DEBUG
			// Output recv_buff
			if (OUTPUT_RANK == rankid) {
				cout << "recv_buff: " << endl;
				for (int i = 0; i < dis_num_item_recvd; ++i) {
					if (recv_buff[i] == -1) {
						cout << endl;
					} else {
						cout << recv_buff[i] << " ";
					}
				}
			}
#endif

			/* Convert recv_buff to matrix vector<unordered_set<T> > disCommunities and copy to send_buff */
			if (message_send_num < (numprocs - 1) && recv_tag != send_dest) {
				dis_num_item_sent = dis_num_item_recvd;
				send_buff = (T *) malloc(dis_num_item_sent * sizeof(T));
			}
			bufCommunities.clear();
			unordered_set<T> community;
			for (int i = 0; i < dis_num_item_recvd; ++i) {
				if (recv_buff[i] == -1) {
					bufCommunities.push_back(community);
					community.clear();
				} else {
					community.insert(recv_buff[i]);
				}

				if (message_send_num < (numprocs - 1)
						&& recv_tag != send_dest) {
					send_buff[i] = recv_buff[i];
				}
			}

			free(recv_buff);

#ifdef MPI_DEBUG
			// Output received communities
			if (OUTPUT_RANK == rankid) {
				int disComSize = bufCommunities.size();
				cout << disComSize << " discovered communities: " << endl;
				for (int i = 0; i < disComSize; ++i) {
					unordered_set<int> community = bufCommunities[i];
					typename unordered_set<int>::const_iterator comIter;
					for (comIter = community.begin();
							comIter != community.end(); ++comIter) {
						cout << *comIter << " ";
					}
					cout << endl;
				}
			}
#endif

			// Compute the FMeasure and NVD from other ranks
			startCompCycle = rdtsc();
			disComNum = bufCommunities.size();
			for (int i = 0; i < realComNum; ++i) {
				unordered_set<T> realCommunity = realCommunities[i];
				double ni = realCommunity.size();
				double maxCommon = -1;
				double maxFMeasure = -1;

				for (int j = 0; j < disComNum; ++j) {
					unordered_set<T> disCommunity = bufCommunities[j];
					double nj = disCommunity.size();
					double nij = this->getIntersectionNum(realCommunity,
							disCommunity);

					if (nij > maxCommon) {
						maxCommon = nij;
					}

					double tmpFMeasure = (2 * nij) / (ni + nj);
					if (tmpFMeasure > maxFMeasure) {
						maxFMeasure = tmpFMeasure;
					}
				}

				if (maxFMeasure > realMaxFMeasures[i]) {
					realMaxFMeasures[i] = maxFMeasure;
				}

				if (maxCommon > realMaxCommons[i]) {
					realMaxCommons[i] = maxCommon;
				}

//			rankNVD += maxCommon;
//			rankFMeasure += ni * maxFMeasure;
			}
			endCompCycle = rdtsc();
			clusteringCycles[0] += endCompCycle - startCompCycle;
		} /* while */
	} // if numprocs>1

// Calculate FMeasure and the first part of NVD
	startCompCycle = rdtsc();
	for (int i = 0; i < realComNum; ++i) {
		double ni = realCommunities[i].size();
		rankFMeasure += ni * realMaxFMeasures[i];
		rankNVD += realMaxCommons[i];
	}
	endCompCycle = rdtsc();
	clusteringCycles[0] += endCompCycle - startCompCycle;

#ifdef MPI_DEBUG
// (2) To get the second part of NVD
	if (OUTPUT_RANK == rankid) {
		int disComSize = disCommunities.size();
		cout << disComSize << " original discovered communities: " << endl;
		for (int i = 0; i < disComSize; ++i) {
			unordered_set<int> community = disCommunities[i];
			typename unordered_set<int>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				cout << *comIter << " ";
			}
			cout << endl;
		}
	}
#endif

	disComNum = disCommunities.size();
// To record the max of each detected community of this rank
	startCompCycle = rdtsc();
	// double disMaxCommons[disComNum];
	double *disMaxCommons = (double *) malloc(disComNum * sizeof(double));
	memset(disMaxCommons, 0, disComNum * sizeof(double));
	for (int i = 0; i < disComNum; ++i) {
		unordered_set<T> disCommunity = disCommunities[i];
		double maxCommon = -1;

		for (int j = 0; j < realComNum; ++j) {
			unordered_set<T> realCommunity = realCommunities[j];
			double nij = this->getIntersectionNum(disCommunity, realCommunity);

			if (nij > maxCommon) {
				maxCommon = nij;
			}
		}

		if (maxCommon > disMaxCommons[i]) {
			disMaxCommons[i] = maxCommon;
		}

		//rankNVD += maxCommon;
	}
	endCompCycle = rdtsc();
	clusteringCycles[0] += endCompCycle - startCompCycle;

// Only if numprocs>1, shall we consider the message passing part
	if (numprocs > 1) {
		message_send_num = 0;

		/* Message send variables */
		send_request = MPI_REQUEST_NULL;
		T *real_send_buff = (T *) malloc(real_num_item_sent * sizeof(T));

		// Initialize the real_send_buff to be the real communities of this rank
		int step = 0;
		for (int i = 0; i < realComNum; ++i) {
			unordered_set<T> community = realCommunities[i];

			typename unordered_set<T>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				real_send_buff[step++] = *comIter;
			}

			real_send_buff[step++] = -1;
		}

#ifdef MPI_DEBUG
		// Output self send_buff
		if (rankid == OUTPUT_RANK) {
			cout << "self real_send_buff: " << endl;
			for (int i = 0; i < real_num_item_sent; ++i) {
				if (real_send_buff[i] == -1) {
					cout << endl;
				} else {
					cout << real_send_buff[i] << " ";
				}
			}
		}
#endif

		/* Message recv variables */
		int real_num_item_recvd = 0;
		T *real_recv_buff;

		/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
		 * next MPI rank at the first step and also the message tag is its rankid.
		 */
		recv_tag = rankid;

		/* When receiving (numprocs - 1) messages, exit */
		while (message_send_num < (numprocs - 1)) {
			startMsgCycle = rdtsc();
			/* At first step, send its own slice of matrix B */
			/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
			if (recv_tag != send_dest) {
				MPI_Isend(real_send_buff, real_num_item_sent, MPI_INT,
						send_dest, recv_tag, MPI_COMM_WORLD, &send_request);
			}

			recv_flag = 0;
			while (!recv_flag) {
				MPI_Iprobe(recv_source, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_flag,
						&status);

				// When ready to receive
				if (recv_flag) {
					// Get the number of item that should be received
					MPI_Get_count(&status, MPI_INT, &real_num_item_recvd);
					// Post a recv to receive the message
					real_recv_buff = (T *) malloc(
							real_num_item_recvd * sizeof(T));
					MPI_Recv(real_recv_buff, real_num_item_recvd, MPI_INT,
							recv_source, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					/* recv_tag is the taskid that owns this slice of matrix B */
					recv_tag = status.MPI_TAG;
				}
			}

			/* wait until the MPI_IRecv and MPI_ISend succeed */
			send_flag = 0;
			while (!send_flag) {
				/* Only when send_request is not null, shall we test its status */
				if (send_request != MPI_REQUEST_NULL) {
					MPI_Test(&send_request, &send_flag, &status);

					/* When finishing sending message, post a recv request */
					if (send_flag) {
						++message_send_num;
						free(real_send_buff);

//						cout << "Second: message_send_num = "
//								<< message_send_num << endl;
					}
				} /* if request */
			} /* while flag */
			endMsgCycle = rdtsc();
			clusteringCycles[1] += endMsgCycle - startMsgCycle;

#ifdef MPI_DEBUG
			// Output recv_buff
			if (OUTPUT_RANK == rankid) {
				cout << "recv_buff: " << endl;
				for (int i = 0; i < real_num_item_recvd; ++i) {
					if (real_recv_buff[i] == -1) {
						cout << endl;
					} else {
						cout << real_recv_buff[i] << " ";
					}
				}
			}
#endif

			/* Convert recv_buff to matrix vector<unordered_set<T> > disCommunities and copy to send_buff */
			if (message_send_num < (numprocs - 1) && recv_tag != send_dest) {
				real_num_item_sent = real_num_item_recvd;
				real_send_buff = (T *) malloc(real_num_item_sent * sizeof(T));
			}
			bufCommunities.clear();
			unordered_set<T> community;
			for (int i = 0; i < real_num_item_recvd; ++i) {
				if (real_recv_buff[i] == -1) {
					bufCommunities.push_back(community);
					community.clear();
				} else {
					community.insert(real_recv_buff[i]);
				}

				if (message_send_num < (numprocs - 1)
						&& recv_tag != send_dest) {
					real_send_buff[i] = real_recv_buff[i];
				}
			}

			free(real_recv_buff);

			// Compute the VI and NMI from other ranks
			startCompCycle = rdtsc();
			realComNum = bufCommunities.size();
			for (int i = 0; i < disComNum; ++i) {
				unordered_set<T> disCommunity = disCommunities[i];
				double maxCommon = -1;

				for (int j = 0; j < realComNum; ++j) {
					unordered_set<T> realCommunity = bufCommunities[j];
					double nij = this->getIntersectionNum(disCommunity,
							realCommunity);

					if (nij > maxCommon) {
						maxCommon = nij;
					}
				}

				if (maxCommon > disMaxCommons[i]) {
					disMaxCommons[i] = maxCommon;
				}
				// rankNVD += maxCommon;
			}
			endCompCycle = rdtsc();
			clusteringCycles[0] += endCompCycle - startCompCycle;

		} /* while */
	}

// Calculate the second part of NVD of this rank
	startCompCycle = rdtsc();
	for (int i = 0; i < disComNum; ++i) {
		rankNVD += disMaxCommons[i];
	}
	endCompCycle = rdtsc();
	clusteringCycles[0] += endCompCycle - startCompCycle;

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankNVD = " << rankNVD << ", rankFMeasure = " << rankFMeasure
		<< endl;
	}
#endif

	double NVD = 0;
	double FMeasure = 0;
	MPI_Allreduce(&rankNVD, &NVD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankFMeasure, &FMeasure, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	FMeasure = FMeasure / numNodes;
	NVD = 1 - NVD / (2 * numNodes);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "FMeasure = " << FMeasure << ", NVD = " << NVD << endl;
	}
#endif

	metrics[0] = FMeasure;
	metrics[1] = NVD;
}
//------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------//
/**
 * Compute F-measure and Normalized Van Dongen metric (NVD): Disjoint
 * community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realCommunities: the real communities of this rank
 * @param disCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (VI and NMI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeClusterMatchingMetric_backup(int rankid, int numprocs,
		const vector<unordered_set<T> >& realCommunities,
		const vector<unordered_set<T> >& disCommunities, double *metrics) {
	int realComNum = realCommunities.size();
	int disComNum = disCommunities.size();

// Get the number of nodes in this rank and also the max number of nodes among real communities
	double rankNumNodes = 0;
	int realRankMaxComSize = -1;
	for (int i = 0; i < realComNum; ++i) {
		int ni = realCommunities[i].size();
		if (ni > realRankMaxComSize) {
			realRankMaxComSize = ni;
		}
		rankNumNodes += ni;
	}

// Get the max number of nodes of a community among all the detected communities of this rank
	int disRankMaxComSize = -1;
	for (int j = 0; j < disComNum; ++j) {
		int nj = disCommunities[j].size();
		if (nj > disRankMaxComSize) {
			disRankMaxComSize = nj;
		}
	}

// Get the total number of nodes of the network
	double numNodes = 0;
	MPI_Allreduce(&rankNumNodes, &numNodes, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

// Get the max number of nodes of a community among all the real communities
	int realMaxComSize = 0;
	MPI_Allreduce(&realRankMaxComSize, &realMaxComSize, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

// Get the max number of real communities that a rank may have
	int realMaxComNum = 0;
	MPI_Allreduce(&realComNum, &realMaxComNum, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

// Get the max number of nodes of a community among all the detected communities
	int disMaxComSize = 0;
	MPI_Allreduce(&disRankMaxComSize, &disMaxComSize, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

// Get the max number of detected communities that a rank may have
	int disMaxComNum = 0;
	MPI_Allreduce(&disComNum, &disMaxComNum, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "numNodes = " << numNodes << ", realMaxComNum = "
		<< realMaxComNum << ", realMaxComSize = " << realMaxComSize
		<< ", disMaxComNum = " << disMaxComNum << ", disMaxComSize = "
		<< disMaxComSize << endl;
	}
#endif

// Compute NVD and FMeasure for this rank.
	double rankNVD = 0;
	double rankFMeasure = 0;

// (1) To get FMeasure and the first part of NVD
// To record the max of each real community of this rank
	double realMaxFMeasures[realComNum];
	double realMaxCommons[realComNum];
	memset(realMaxFMeasures, 0, realComNum * sizeof(double));
	memset(realMaxCommons, 0, realComNum * sizeof(double));
	for (int i = 0; i < realComNum; ++i) {
		unordered_set<T> realCommunity = realCommunities[i];
		double ni = realCommunity.size();
		double maxCommon = -1;
		double maxFMeasure = -1;

		for (int j = 0; j < disComNum; ++j) {
			unordered_set<T> disCommunity = disCommunities[j];
			double nj = disCommunity.size();
			double nij = this->getIntersectionNum(realCommunity, disCommunity);

			if (nij > maxCommon) {
				maxCommon = nij;
			}

			double tmpFMeasure = (2 * nij) / (ni + nj);
			if (tmpFMeasure > maxFMeasure) {
				maxFMeasure = tmpFMeasure;
			}
		}

		if (maxFMeasure > realMaxFMeasures[i]) {
			realMaxFMeasures[i] = maxFMeasure;
		}

		if (maxCommon > realMaxCommons[i]) {
			realMaxCommons[i] = maxCommon;
		}
//		rankNVD += maxCommon;
//		rankFMeasure += ni * maxFMeasure;
	}

	int message_send_num = 0;
	int send_flag, recv_flag = 0;
	MPI_Status status;

// The total number of elements will be sent or received on detected communities
	unsigned int count = disMaxComNum * disMaxComSize;

	/* Message send and recv variables */
	MPI_Request recv_request = MPI_REQUEST_NULL;
	int recv_source = (rankid + 1) % numprocs;

	/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
	 * next MPI rank at the first step and also the message tag is its rankid.
	 */
	int recv_tag = rankid;

// Save the received communities
	vector<unordered_set<T> > bufCommunities;

	/* Message send variables */
	MPI_Request send_request = MPI_REQUEST_NULL;
	int send_dest = (rankid + numprocs - 1) % numprocs;

// Only if numprocs>1, shall we consider the message passing part
	if (numprocs > 1) {
		//T **recv_buff = malloc2d(disMaxComNum, disMaxComSize);
		T recv_buff[disMaxComNum][disMaxComSize];

		//T **send_buff = malloc2d(disMaxComNum, disMaxComSize);
		T send_buff[disMaxComNum][disMaxComSize];

		// Initialize the send_buff to be the detected communities of this rank
		for (int i = 0; i < disComNum; ++i) {
			unordered_set<T> community = disCommunities[i];
			int step = 0;
			typename unordered_set<T>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				send_buff[i][step] = *comIter;
				++step;
			}

			// Fill the remaining elements with END_OF_COMMUNITY
			for (int j = step; j < disMaxComSize; ++j) {
				send_buff[i][j] = END_OF_COMMUNITY;
			}
		}

		// Fill the remaining rows with END_OF_COMMUNITY
		for (int i = disComNum; i < disMaxComNum; ++i) {
			for (int j = 0; j < disMaxComSize; ++j) {
				send_buff[i][j] = END_OF_COMMUNITY;
			}
		}

#ifdef MPI_DEBUG
		// Output self send_buff
		if (rankid == OUTPUT_RANK) {
			cout << "self send_buff: " << endl;
			for (int i = 0; i < disMaxComNum; ++i) {
				for (int j = 0; j < disMaxComSize; ++j) {
					cout << send_buff[i][j] << " ";
				}
				cout << endl;
			}
		}
#endif

		/* When receiving (numprocs - 1) messages, exit */
		while (message_send_num < (numprocs - 1)) {
			/* First, post a recv request to prevent deadlock */
			MPI_Irecv(&recv_buff[0][0], count, MPI_INT, recv_source,
					MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);

			/* At first step, send its own slice of matrix B */
			/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
			if (recv_tag != send_dest) {
				MPI_Isend(&send_buff[0][0], count, MPI_INT, send_dest, recv_tag,
						MPI_COMM_WORLD, &send_request);
			}

			/* wait until the MPI_IRecv and MPI_ISend succeed */
			send_flag = 0;
			recv_flag = 0;
			while (!recv_flag || !send_flag) {
				/* Only when recv_request is not null, shall we test its status */
				if (recv_request != MPI_REQUEST_NULL) {
					MPI_Test(&recv_request, &recv_flag, &status);

					/* When finishing recving message, post a send request */
					if (recv_flag) {
						/* recv_tag is the taskid that owns this slice of matrix B */
						recv_tag = status.MPI_TAG;
					} /* if flag */
				} /* if request*/

				/* Only when send_request is not null, shall we test its status */
				if (send_request != MPI_REQUEST_NULL) {
					MPI_Test(&send_request, &send_flag, &status);

					/* When finishing sending message, post a recv request */
					if (send_flag) {
						++message_send_num;
					}
				} /* if request */
			} /* while flag */

#ifdef MPI_DEBUG
			// Output recv_buff
			if (OUTPUT_RANK == rankid) {
				cout << "recv_buff: " << endl;
				for (int i = 0; i < disMaxComNum; ++i) {
					for (int j = 0; j < disMaxComSize; ++j) {
						cout << recv_buff[i][j] << " ";
					}
					cout << endl;
				}
			}
#endif

			// Convert recv_buff to vector<unordered_set<T> > disCommunities and copy to send_buff
			bufCommunities.clear();
			for (int i = 0; i < disMaxComNum; ++i) {
				// If this line does not correspond to a community, ignore it
				if (recv_buff[i][0] == END_OF_COMMUNITY) {
					for (int j = 0; j < disMaxComSize; ++j) {
						send_buff[i][j] = recv_buff[i][j];
					}
					continue;
				}
				unordered_set<T> community;
				for (int j = 0; j < disMaxComSize; ++j) {
					int nodeId = recv_buff[i][j];
					send_buff[i][j] = nodeId;
					if (nodeId != END_OF_COMMUNITY) {
						community.insert(nodeId);
					}
				}

				bufCommunities.push_back(community);
			}

#ifdef MPI_DEBUG
			// Output received communities
			if (OUTPUT_RANK == rankid) {
				int disComNum = bufCommunities.size();
				cout << disComNum << " discovered communities: " << endl;
				for (int i = 0; i < disComNum; ++i) {
					unordered_set<int> community = bufCommunities[i];
					typename unordered_set<int>::const_iterator comIter;
					for (comIter = community.begin();
							comIter != community.end(); ++comIter) {
						cout << *comIter << " ";
					}
					cout << endl;
				}
			}
#endif

			// Compute the VI and NMI from other ranks
			disComNum = bufCommunities.size();
			for (int i = 0; i < realComNum; ++i) {
				unordered_set<T> realCommunity = realCommunities[i];
				double ni = realCommunity.size();
				double maxCommon = -1;
				double maxFMeasure = -1;

				for (int j = 0; j < disComNum; ++j) {
					unordered_set<T> disCommunity = bufCommunities[j];
					double nj = disCommunity.size();
					double nij = this->getIntersectionNum(realCommunity,
							disCommunity);

					if (nij > maxCommon) {
						maxCommon = nij;
					}

					double tmpFMeasure = (2 * nij) / (ni + nj);
					if (tmpFMeasure > maxFMeasure) {
						maxFMeasure = tmpFMeasure;
					}
				}

				if (maxFMeasure > realMaxFMeasures[i]) {
					realMaxFMeasures[i] = maxFMeasure;
				}

				if (maxCommon > realMaxCommons[i]) {
					realMaxCommons[i] = maxCommon;
				}

//			rankNVD += maxCommon;
//			rankFMeasure += ni * maxFMeasure;
			}
		} /* while */
	} // if numprocs>1

// Calculate FMeasure and the first part of NVD
	for (int i = 0; i < realComNum; ++i) {
		double ni = realCommunities[i].size();
		rankFMeasure += ni * realMaxFMeasures[i];
		rankNVD += realMaxCommons[i];
	}

#ifdef MPI_DEBUG
// (2) To get the second part of NVD
	if (OUTPUT_RANK == rankid) {
		int disComSize = disCommunities.size();
		cout << disComSize << " original discovered communities: " << endl;
		for (int i = 0; i < disComSize; ++i) {
			unordered_set<int> community = disCommunities[i];
			typename unordered_set<int>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				cout << *comIter << " ";
			}
			cout << endl;
		}
	}
#endif

	disComNum = disCommunities.size();
// To record the max of each detected community of this rank
	double disMaxCommons[disComNum];
	memset(disMaxCommons, 0, disComNum * sizeof(double));
	for (int i = 0; i < disComNum; ++i) {
		unordered_set<T> disCommunity = disCommunities[i];
		double maxCommon = -1;

		for (int j = 0; j < realComNum; ++j) {
			unordered_set<T> realCommunity = realCommunities[j];
			double nij = this->getIntersectionNum(disCommunity, realCommunity);

			if (nij > maxCommon) {
				maxCommon = nij;
			}
		}

		if (maxCommon > disMaxCommons[i]) {
			disMaxCommons[i] = maxCommon;
		}

		//rankNVD += maxCommon;
	}

// Only if numprocs>1, shall we consider the message passing part
	if (numprocs > 1) {
		message_send_num = 0;
		count = realMaxComNum * realMaxComSize;

		/* Message send variables */
		send_request = MPI_REQUEST_NULL;
		//T **real_send_buff = malloc2d(realMaxComNum, maxComSize);
		T real_send_buff[realMaxComNum][realMaxComSize];

		// Initialize the real_send_buff to be the real communities of this rank
		for (int i = 0; i < realComNum; ++i) {
			unordered_set<T> community = realCommunities[i];
			int step = 0;
			typename unordered_set<T>::const_iterator comIter;
			for (comIter = community.begin(); comIter != community.end();
					++comIter) {
				real_send_buff[i][step] = *comIter;
				++step;
			}

			// Fill the remaining elements with END_OF_COMMUNITY
			for (int j = step; j < realMaxComSize; ++j) {
				real_send_buff[i][j] = END_OF_COMMUNITY;
			}
		}

		// Fill the remaining rows with END_OF_COMMUNITY
		for (int i = realComNum; i < realMaxComNum; ++i) {
			for (int j = 0; j < realMaxComSize; ++j) {
				real_send_buff[i][j] = END_OF_COMMUNITY;
			}
		}

#ifdef MPI_DEBUG
		// Output self send_buff
		if (rankid == OUTPUT_RANK) {
			cout << "self real_send_buff: " << endl;
			for (int i = 0; i < realMaxComNum; ++i) {
				for (int j = 0; j < realMaxComSize; ++j) {
					cout << real_send_buff[i][j] << " ";
				}
				cout << endl;
			}
		}
#endif

		/* Message recv variables */
		recv_request = MPI_REQUEST_NULL;
		//T **real_recv_buff = malloc2d(realMaxComNum, realMaxComSize);
		T real_recv_buff[realMaxComNum][realMaxComSize];

		/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
		 * next MPI rank at the first step and also the message tag is its rankid.
		 */
		recv_tag = rankid;

		/* When receiving (numprocs - 1) messages, exit */
		while (message_send_num < (numprocs - 1)) {
			/* First, post a recv request to prevent deadlock */
			MPI_Irecv(&real_recv_buff[0][0], count, MPI_INT, recv_source,
					MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);

			/* At first step, send its own slice of matrix B */
			/* At other steps, don't send the slice of B to the MPI rank that owns this part of B */
			if (recv_tag != send_dest) {
				MPI_Isend(&real_send_buff[0][0], count, MPI_INT, send_dest,
						recv_tag, MPI_COMM_WORLD, &send_request);
			}

			/* wait until the MPI_IRecv and MPI_ISend succeed */
			send_flag = 0;
			recv_flag = 0;
			while (!recv_flag || !send_flag) {
				/* Only when recv_request is not null, shall we test its status */
				if (recv_request != MPI_REQUEST_NULL) {
					MPI_Test(&recv_request, &recv_flag, &status);

					/* When finishing recving message, post a send request */
					if (recv_flag) {
						/* recv_tag is the taskid that owns this slice of matrix B */
						recv_tag = status.MPI_TAG;
					} /* if flag */
				} /* if request*/

				/* Only when send_request is not null, shall we test its status */
				if (send_request != MPI_REQUEST_NULL) {
					MPI_Test(&send_request, &send_flag, &status);

					/* When finishing sending message, post a recv request */
					if (send_flag) {
						++message_send_num;
					}
				} /* if request */
			} /* while flag */

#ifdef MPI_DEBUG
			// Output recv_buff
			if (OUTPUT_RANK == rankid) {
				cout << "real_recv_buff: " << endl;
				for (int i = 0; i < realMaxComNum; ++i) {
					for (int j = 0; j < realMaxComSize; ++j) {
						cout << real_recv_buff[i][j] << " ";
					}
					cout << endl;
				}
			}
#endif

			/* Convert recv_buff to matrix vector<unordered_set<T> > disCommunities */
			bufCommunities.clear();
			for (int i = 0; i < realMaxComNum; ++i) {
				// If this line does not correspond to a community, ignore it
				if (real_recv_buff[i][0] == END_OF_COMMUNITY) {
					for (int j = 0; j < realMaxComSize; ++j) {
						real_send_buff[i][j] = real_recv_buff[i][j];
					}
					continue;
				}
				unordered_set<T> community;
				for (int j = 0; j < realMaxComSize; ++j) {
					int nodeId = real_recv_buff[i][j];
					real_send_buff[i][j] = nodeId;
					if (nodeId != END_OF_COMMUNITY) {
						community.insert(nodeId);
					}
				}

				bufCommunities.push_back(community);
			}

#ifdef MPI_DEBUG
			// Output received communities
			if (OUTPUT_RANK == rankid) {
				int disComNum = bufCommunities.size();
				cout << disComNum << " real communities: " << endl;
				for (int i = 0; i < disComNum; ++i) {
					unordered_set<int> community = bufCommunities[i];
					typename unordered_set<int>::const_iterator comIter;
					for (comIter = community.begin();
							comIter != community.end(); ++comIter) {
						cout << *comIter << " ";
					}
					cout << endl;
				}
			}
#endif

			// Compute the VI and NMI from other ranks
			realComNum = bufCommunities.size();
			for (int i = 0; i < disComNum; ++i) {
				unordered_set<T> disCommunity = disCommunities[i];
				double maxCommon = -1;

				for (int j = 0; j < realComNum; ++j) {
					unordered_set<T> realCommunity = bufCommunities[j];
					double nij = this->getIntersectionNum(disCommunity,
							realCommunity);

					if (nij > maxCommon) {
						maxCommon = nij;
					}
				}

				if (maxCommon > disMaxCommons[i]) {
					disMaxCommons[i] = maxCommon;
				}
				// rankNVD += maxCommon;
			}
		} /* while */
	}

// Calculate the second part of NVD of this rank
	for (int i = 0; i < disComNum; ++i) {
		rankNVD += disMaxCommons[i];
	}

// free memory
//free2d((T **) real_send_buff);
//free2d((T **) real_recv_buff);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankNVD = " << rankNVD << ", rankFMeasure = " << rankFMeasure
		<< endl;
	}
#endif

	double NVD = 0;
	double FMeasure = 0;
	MPI_Allreduce(&rankNVD, &NVD, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankFMeasure, &FMeasure, 1, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);
	FMeasure = FMeasure / numNodes;
	NVD = 1 - NVD / (2 * numNodes);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "FMeasure = " << FMeasure << ", NVD = " << NVD << endl;
	}
#endif

	metrics[0] = FMeasure;
	metrics[1] = NVD;
}
//------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------//
/**
 * Compute Rand Index, Adjusted Rand Index, and Jaccard Index: Disjoint
 * community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realMapCommunities: the real communities of this rank
 * @param disMapCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (RI, ARI, and JI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeIndexMetric(int rankid, int numprocs,
		unordered_map<T, T>& realMapCommunities,
		unordered_map<T, T>& disMapCommunities, double *metrics,
		double *indexCycles) {
	double rankA11 = 0;
	double rankA00 = 0;
	double rankA10 = 0;
	double rankA01 = 0;
	int rankNodeNum = realMapCommunities.size();

	// Compute RI, ARI, JI of this rank
	double startCompCycle = rdtsc();
	typename unordered_map<T, T>::const_iterator realIter;
	typename unordered_map<T, T>::const_iterator iIter;
	typename unordered_map<T, T>::const_iterator jIter;
	for (iIter = realMapCommunities.begin(); iIter != realMapCommunities.end();
			++iIter) {
		int inode = iIter->first;
		for (jIter = realMapCommunities.begin();
				jIter != realMapCommunities.end(); ++jIter) {
			int jnode = jIter->first;
			if (inode != jnode) {
				int icom = disMapCommunities[inode];
				int jcom = disMapCommunities[jnode];
				bool dflag = true;
				if (icom != jcom) {
					dflag = false;
				}
				icom = realMapCommunities[inode];
				jcom = realMapCommunities[jnode];
				bool rflag = true;
				if (icom != jcom) {
					rflag = false;
				}

				if (dflag && rflag) {
					rankA11 += 1;
				} else if (rflag && !dflag) {
					rankA10 += 1;
				} else if (!rflag && dflag) {
					rankA01 += 1;
				} else if (!rflag && !dflag) {
					rankA00 += 1;
				}
			}
		}
	} // for
	double endCompCycle = rdtsc();
	indexCycles[0] += endCompCycle - startCompCycle;

// Only if numprocs>1, shall we go to the message passing part
	if (numprocs > 1) {
		double startMsgCycle, endMsgCycle;
		// Message send variables
		unsigned int message_send_num = 0;
		unsigned int num_item_sent = 2 * rankNodeNum;
		int send_flag = 0, dis_send_flag = 0;
		MPI_Request send_request = MPI_REQUEST_NULL;
		MPI_Request dis_send_request = MPI_REQUEST_NULL;
		int send_dest = (rankid + numprocs - 1) % numprocs;
		/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
		 * next MPI rank at the first step and also the message tag is its rankid.
		 */
		int send_tag = rankid;

		// Detected nodes buff with community information
		int index = 0;

#ifdef MPI_DEBUG
		if (OUTPUT_RANK == rankid) {
			cout << "BEGIN> " << rankNodeNum << " detected nodes: " << endl;
		}
#endif

		T *dis_send_buff = (T *) malloc(num_item_sent * sizeof(T));
		typename unordered_map<T, T>::const_iterator disIter;
		for (disIter = disMapCommunities.begin();
				disIter != disMapCommunities.end(); ++disIter) {
			dis_send_buff[index++] = disIter->first;
			dis_send_buff[index++] = disIter->second;

#ifdef MPI_DEBUG
			if (OUTPUT_RANK == rankid) {
				cout << dis_send_buff[index - 2] << "\t"
				<< dis_send_buff[index - 1] << endl;
			}
#endif
		}

		/* Message recv variables */
		unsigned int num_item_recvd = 0;
		int recv_flag = 0, dis_recv_flag = 0;
		MPI_Request recv_request = MPI_REQUEST_NULL;
		MPI_Request dis_recv_request = MPI_REQUEST_NULL;
		int recv_source = (rankid + 1) % numprocs;
		int recv_tag;
		// Save the received node communities
		unordered_map<T, T> bufMapCommunities;

		MPI_Status status;

		/* When requirements meet, exit */
		while (message_send_num < (numprocs - 1)) {
#ifdef MPI_DEBUG
			cout << "while: rankid = " << rankid << ", message_send_num = "
			<< message_send_num << endl;
#endif
			startMsgCycle = rdtsc();
			// First, post a recv request to get the number of nodes will be received
			MPI_Irecv(&num_item_recvd, 1, MPI_INT, recv_source, MPI_ANY_TAG,
					MPI_COMM_WORLD, &recv_request);

			if (send_tag != send_dest) {
				// Then, post a send request to send the number of nodes will be sent
				MPI_Isend(&num_item_sent, 1, MPI_INT, send_dest, send_tag,
						MPI_COMM_WORLD, &send_request);
			}

			// Test recv_request
			recv_flag = 0;
			while (!recv_flag) {
				/* Only when recv_request is not null, shall we test its status */
				if (recv_request != MPI_REQUEST_NULL) {
					MPI_Test(&recv_request, &recv_flag, &status);

					/* When finishing recving message, post a send request */
					if (recv_flag) {
						/* recv_tag is the taskid that owns this slice of matrix B */
						recv_tag = status.MPI_TAG;
#ifdef MPI_DEBUG
						cout << "VAR: rankid = " << rankid << ", recv_tag = "
						<< recv_tag << endl;
#endif
					} /* if flag */
				} /* if request*/
			}

			// Test send_request
			if (send_tag != send_dest) {
				send_flag = 0;
				while (!send_flag) {
					/* Only when recv_request is not null, shall we test its status */
					if (send_request != MPI_REQUEST_NULL) {
						MPI_Test(&send_request, &send_flag, &status);
					} /* if request*/
				}

#ifdef MPI_DEBUG
				if (OUTPUT_RANK == rankid) {
					cout << "rankid = " << rankid << ", num_item_sent = "
					<< num_item_sent << endl;
				}
#endif
			}

			// Receive detected community nodes
			/* Second, post a recv request to get detected community nodes */
			T *dis_recv_buff = (T *) malloc(num_item_recvd * sizeof(T));
			MPI_Irecv(dis_recv_buff, num_item_recvd, MPI_INT, recv_source,
					MPI_ANY_TAG, MPI_COMM_WORLD, &dis_recv_request);

			/* At first step, send its own community nodes */
			/* At other steps, don't send community nodes to the MPI rank that owns them */
			if (send_tag != send_dest) {
				MPI_Isend(dis_send_buff, num_item_sent, MPI_INT, send_dest,
						send_tag, MPI_COMM_WORLD, &dis_send_request);
			}

			// Test receive request
			dis_recv_flag = 0;
			while (!dis_recv_flag) {
				/* Only when send_request is not null, shall we test its status */
				if (dis_recv_request != MPI_REQUEST_NULL) {
					MPI_Test(&dis_recv_request, &dis_recv_flag, &status);

					/* When finishing sending message, post a recv request */
					if (dis_recv_flag) {
						recv_tag = status.MPI_TAG;
#ifdef MPI_DEBUG
						cout << "BUFF: rankid = " << rankid << ", recv_tag = "
						<< recv_tag << endl;
#endif
					}
				} /* if request */
			} /* while flag */

#ifdef MPI_DEBUG
			cout << "rankid = " << rankid << " succeed in receiving message!"
			<< endl;
#endif

			// Test send request
			if (send_tag != send_dest) {
				/* wait until MPI_ISend succeed */
				dis_send_flag = 0;
				while (!dis_send_flag) {
					/* Only when send_request is not null, shall we test its status */
					if (dis_send_request != MPI_REQUEST_NULL) {
						MPI_Test(&dis_send_request, &dis_send_flag, &status);
					} /* if request */
				} /* while flag */

				++message_send_num;

#ifdef MPI_DEBUG
				cout << "rankid = " << rankid
				<< " succeed in sending message!, message_send_num = "
				<< message_send_num << endl;
#endif

				free(dis_send_buff);
			}
			endMsgCycle = rdtsc();
			indexCycles[1] += endMsgCycle - startMsgCycle;

			// Construct node map communities from received array
#ifdef MPI_DEBUG
			if (OUTPUT_RANK == rankid) {
				cout << num_item_recvd / 2 << " received nodes: " << endl;
			}
#endif
			send_tag = recv_tag;
			bufMapCommunities.clear();
			// If should send message next turn, copy recv_buffs to send_buffs
			if (message_send_num < (numprocs - 1) && send_tag != send_dest) {
				num_item_sent = num_item_recvd;
				dis_send_buff = (T *) malloc(num_item_sent * sizeof(T));
			}
			for (int i = 0; i < num_item_recvd; i += 2) {
				bufMapCommunities.insert(
						make_pair(dis_recv_buff[i], dis_recv_buff[i + 1]));
				if (message_send_num < (numprocs - 1)
						&& send_tag != send_dest) {
					dis_send_buff[i] = dis_recv_buff[i];
					dis_send_buff[i + 1] = dis_recv_buff[i + 1];
				}
#ifdef MPI_DEBUG
				if (OUTPUT_RANK == rankid) {
					cout << dis_recv_buff[i] << "\t" << dis_recv_buff[i + 1]
					<< endl;
				}
#endif
			}

			// Compute RI, ARI, JI from other ranks
			startCompCycle = rdtsc();
			for (realIter = realMapCommunities.begin();
					realIter != realMapCommunities.end(); ++realIter) {
				int inode = realIter->first;
				for (int j = 0; j < num_item_recvd; j += 2) {
					int jnode = dis_recv_buff[j];
					int icom = disMapCommunities[inode];
					int jcom = bufMapCommunities[jnode];
					if (icom == jcom) {
						rankA01 += 1;
					} else {
						rankA00 += 1;
					}
				}
			} // for
			endCompCycle = rdtsc();
			indexCycles[0] += endCompCycle - startCompCycle;

			free(dis_recv_buff);

#ifdef MPI_DEBUG
			cout << "rankid = " << rankid << " finishes computation!" << endl;
#endif
		}
		/* while */
	} // if numprocs>1

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankid = " << rankid << ", rankA11 = " << rankA11
		<< ", rankA00 = " << rankA00 << ", rankA10 = " << rankA10
		<< ", rankA01 = " << rankA01 << endl;
	}
#endif

	double a11 = 0;
	double a00 = 0;
	double a10 = 0;
	double a01 = 0;
	MPI_Allreduce(&rankA11, &a11, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA00, &a00, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA10, &a10, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA01, &a01, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double randIndex = (a11 + a00) / (a11 + a10 + a01 + a00);
	double m = a11 + a10 + a01 + a00;
	double adjustRandIndex = (a11 - ((a11 + a10) * (a11 + a01)) / m)
			/ ((a11 + a10 + a11 + a01) / 2 - ((a11 + a10) * (a11 + a01)) / m);
	double jaccardIndex = a11 / (a11 + a10 + a01);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "RI = " << randIndex << ", ARI = " << adjustRandIndex
		<< ", JI = " << jaccardIndex << endl;
	}
#endif

	metrics[0] = randIndex;
	metrics[1] = adjustRandIndex;
	metrics[2] = jaccardIndex;
}
//--------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------//
/**
 * Compute Rand Index, Adjusted Rand Index, and Jaccard Index: Disjoint
 * community quality
 *
 * @param rankid: the id of this MPI rank
 * @param numprocs: the total number of MPI ranks
 * @param realMapCommunities: the real communities of this rank
 * @param disMapCommunities: the detected communities of this rank
 * @param metrics[]: the calculated metrics (RI, ARI, and JI)
 * @return
 */
template<class T>
void MPIMetric<T>::computeIndexMetric_backup(int rankid, int numprocs,
		unordered_map<T, T>& realMapCommunities,
		unordered_map<T, T>& disMapCommunities, double *metrics,
		double *indexCycles) {
	double rankA11 = 0;
	double rankA00 = 0;
	double rankA10 = 0;
	double rankA01 = 0;
	int rankNodeNum = realMapCommunities.size();

// Real nodes without community information
	double startCompCycle = rdtsc();
	T *realNodes = (T *) malloc(rankNodeNum * sizeof(T));
	typename unordered_map<T, T>::const_iterator realIter;
	int index = 0;
	for (realIter = realMapCommunities.begin();
			realIter != realMapCommunities.end(); ++realIter) {
		realNodes[index++] = realIter->first;
	}

// Compute RI, ARI, JI of this rank
	for (int i = 0; i < rankNodeNum; ++i) {
		int inode = realNodes[i];
		for (int j = i + 1; j < rankNodeNum; ++j) {
			int jnode = realNodes[j];
			if (inode != jnode) {
				int icom = disMapCommunities[inode];
				int jcom = disMapCommunities[jnode];
				bool dflag = true;
				if (icom != jcom) {
					dflag = false;
				}
				icom = realMapCommunities[inode];
				jcom = realMapCommunities[jnode];
				bool rflag = true;
				if (icom != jcom) {
					rflag = false;
				}

				if (dflag && rflag) {
					rankA11 += 1;
				} else if (rflag && !dflag) {
					rankA10 += 1;
				} else if (!rflag && dflag) {
					rankA01 += 1;
				} else if (!rflag && !dflag) {
					rankA00 += 1;
				}
			}
		}
	} // for

	free(realNodes);
	double endCompCycle = rdtsc();
	indexCycles[0] += endCompCycle - startCompCycle;

// Only if numprocs>1, shall we go to the message passing part
	if (numprocs > 1) {
		double startMsgCycle, endMsgCycle;
		unsigned int SEND_MSG_NUM = 0;
		unsigned int REV_MSG_NUM = numprocs - 1;
		if (rankid > 0) {
			SEND_MSG_NUM = numprocs - rankid;
			REV_MSG_NUM = numprocs - 1 - rankid;
		}

		// Message send variables
		unsigned int msg_send_num = 0;
		unsigned int num_item_sent = 2 * rankNodeNum;
		int send_flag = 0, dis_send_flag = 0;
		MPI_Request send_request = MPI_REQUEST_NULL;
		MPI_Request dis_send_request = MPI_REQUEST_NULL;
		int send_dest = (rankid + numprocs - 1) % numprocs;
		/* Set to be rankid to ensure that every MPI rank can send its own community matrix to its
		 * next MPI rank at the first step and also the message tag is its rankid.
		 */
		int send_tag = rankid;

		T *dis_send_buff;
		if (msg_send_num < SEND_MSG_NUM) {
			// Detected nodes buff with community information
			index = 0;

#ifdef MPI_DEBUG
			if (OUTPUT_RANK == rankid) {
				cout << "BEGIN> " << rankNodeNum << " detected nodes: " << endl;
			}
#endif

			dis_send_buff = (T *) malloc(num_item_sent * sizeof(T));
			typename unordered_map<T, T>::const_iterator disIter;
			for (disIter = disMapCommunities.begin();
					disIter != disMapCommunities.end(); ++disIter) {
				dis_send_buff[index++] = disIter->first;
				dis_send_buff[index++] = disIter->second;

#ifdef MPI_DEBUG
				if (OUTPUT_RANK == rankid) {
					cout << dis_send_buff[index - 2] << "\t"
					<< dis_send_buff[index - 1] << endl;
				}
#endif
			}
		}

		/* Message recv variables */
		unsigned int msg_rev_num = 0;
		unsigned int num_item_recvd = 0;
		int recv_flag = 0, dis_recv_flag = 0;
		MPI_Request recv_request = MPI_REQUEST_NULL;
		MPI_Request dis_recv_request = MPI_REQUEST_NULL;
		int recv_source = (rankid + 1) % numprocs;
		int recv_tag;
		// Save the received node communities
		unordered_map<T, T> bufMapCommunities;

		MPI_Status status;

		/* When requirements meet, exit */
		while (msg_send_num < SEND_MSG_NUM || msg_rev_num < REV_MSG_NUM) {
#ifdef MPI_DEBUG
			cout << "while: rankid = " << rankid << ", SEND_MSG_NUM = "
			<< SEND_MSG_NUM << ", msg_send_num = " << msg_send_num
			<< ", REV_MSG_NUM = " << REV_MSG_NUM << ", msg_rev_num = "
			<< msg_rev_num << endl;
#endif
			startMsgCycle = rdtsc();
			if (msg_rev_num < REV_MSG_NUM) {
				// First, post a recv request to get the number of nodes will be received
				MPI_Irecv(&num_item_recvd, 1, MPI_INT, recv_source, MPI_ANY_TAG,
						MPI_COMM_WORLD, &recv_request);
			}
			if (msg_send_num < SEND_MSG_NUM && send_tag != send_dest) {
				// Then, post a send request to send the number of nodes will be sent
				MPI_Isend(&num_item_sent, 1, MPI_INT, send_dest, send_tag,
						MPI_COMM_WORLD, &send_request);
			}

			// Test recv_request
			if (msg_rev_num < REV_MSG_NUM) {
				recv_flag = 0;
				while (!recv_flag) {
					/* Only when recv_request is not null, shall we test its status */
					if (recv_request != MPI_REQUEST_NULL) {
						MPI_Test(&recv_request, &recv_flag, &status);

						/* When finishing recving message, post a send request */
						if (recv_flag) {
							/* recv_tag is the taskid that owns this slice of matrix B */
							recv_tag = status.MPI_TAG;
#ifdef MPI_DEBUG
							cout << "VAR: rankid = " << rankid
							<< ", recv_tag = " << recv_tag << endl;
#endif
						} /* if flag */
					} /* if request*/
				}
			}

			// Test send_request
			if (msg_send_num < SEND_MSG_NUM && send_tag != send_dest) {
				send_flag = 0;
				while (!send_flag) {
					/* Only when recv_request is not null, shall we test its status */
					if (send_request != MPI_REQUEST_NULL) {
						MPI_Test(&send_request, &send_flag, &status);
					} /* if request*/
				}

#ifdef MPI_DEBUG
				if (OUTPUT_RANK == rankid) {
					cout << "rankid = " << rankid << ", num_item_sent = "
					<< num_item_sent << endl;
				}
#endif
			}

			// Receive detected community nodes
			T *dis_recv_buff;
			if (msg_rev_num < REV_MSG_NUM) {
				/* Second, post a recv request to get detected community nodes */
				dis_recv_buff = (T *) malloc(num_item_recvd * sizeof(T));
				MPI_Irecv(dis_recv_buff, num_item_recvd, MPI_INT, recv_source,
						MPI_ANY_TAG, MPI_COMM_WORLD, &dis_recv_request);
			}

			/* At first step, send its own community nodes */
			/* At other steps, don't send community nodes to the MPI rank that owns them */
			if (msg_send_num < SEND_MSG_NUM && send_tag != send_dest) {
				MPI_Isend(dis_send_buff, num_item_sent, MPI_INT, send_dest,
						send_tag, MPI_COMM_WORLD, &dis_send_request);
			}

			// Test receive request
			if (msg_rev_num < REV_MSG_NUM) {
				dis_recv_flag = 0;
				while (!dis_recv_flag) {
					/* Only when send_request is not null, shall we test its status */
					if (dis_recv_request != MPI_REQUEST_NULL) {
						MPI_Test(&dis_recv_request, &dis_recv_flag, &status);

						/* When finishing sending message, post a recv request */
						if (dis_recv_flag) {
							recv_tag = status.MPI_TAG;
#ifdef MPI_DEBUG
							cout << "BUFF: rankid = " << rankid
							<< ", recv_tag = " << recv_tag << endl;
#endif
						}
					} /* if request */
				} /* while flag */

#ifdef MPI_DEBUG
				cout << "rankid = " << rankid
				<< " succeed in receiving message!" << endl;
#endif
			}

			// Test send request
			if (msg_send_num < SEND_MSG_NUM && send_tag != send_dest) {
				/* wait until MPI_ISend succeed */
				dis_send_flag = 0;
				while (!dis_send_flag) {
					/* Only when send_request is not null, shall we test its status */
					if (dis_send_request != MPI_REQUEST_NULL) {
						MPI_Test(&dis_send_request, &dis_send_flag, &status);
					} /* if request */
				} /* while flag */

				++msg_send_num;

#ifdef MPI_DEBUG
				cout << "rankid = " << rankid
				<< " succeed in sending message!, msg_send_num = "
				<< msg_send_num << endl;
#endif

				free(dis_send_buff);
			}
			endMsgCycle = rdtsc();
			indexCycles[1] += endMsgCycle - startMsgCycle;

#ifdef MPI_DEBUG
			if (OUTPUT_RANK == rankid) {
				cout << "rankid = " << rankid << ", SEND_MSG_NUM = "
				<< SEND_MSG_NUM << ", msg_send_num = " << msg_send_num
				<< ", REV_MSG_NUM = " << REV_MSG_NUM
				<< ", msg_rev_num = " << msg_rev_num << endl;
			}
#endif

			if (msg_rev_num < REV_MSG_NUM) {
				// Construct node map communities from received array
#ifdef MPI_DEBUG
				if (OUTPUT_RANK == rankid) {
					cout << num_item_recvd / 2 << " received nodes: " << endl;
				}
#endif
				bufMapCommunities.clear();
				for (int i = 0; i < num_item_recvd; i += 2) {
					bufMapCommunities.insert(
							make_pair(dis_recv_buff[i], dis_recv_buff[i + 1]));
#ifdef MPI_DEBUG
					if (OUTPUT_RANK == rankid) {
						cout << dis_recv_buff[i] << "\t" << dis_recv_buff[i + 1]
						<< endl;
					}
#endif
				}

				// Compute RI, ARI, JI from other ranks
				startCompCycle = rdtsc();
				for (realIter = realMapCommunities.begin();
						realIter != realMapCommunities.end(); ++realIter) {
					int inode = realIter->first;
					for (int j = 0; j < num_item_recvd; j += 2) {
						int jnode = dis_recv_buff[j];
						int icom = disMapCommunities[inode];
						int jcom = bufMapCommunities[jnode];
						if (icom == jcom) {
							rankA01 += 1;
						} else {
							rankA00 += 1;
						}
					}
				} // for
				endCompCycle = rdtsc();
				indexCycles[0] += endCompCycle - startCompCycle;

				send_tag = recv_tag;
				// If should send message next turn, copy recv_buffs to send_buffs
				if (msg_send_num < SEND_MSG_NUM && send_tag != send_dest) {
					num_item_sent = num_item_recvd;
					dis_send_buff = (T *) malloc(num_item_sent * sizeof(T));
					for (int i = 0; i < num_item_sent; ++i) {
						dis_send_buff[i] = dis_recv_buff[i];
					}
				}

				free(dis_recv_buff);
				++msg_rev_num;

#ifdef MPI_DEBUG
				cout << "rankid = " << rankid
				<< " finishes computation!, msg_rev_num = "
				<< msg_rev_num << endl;
#endif
			} // msg_rev_num < REV_MSG_NUM
		}
		/* while */
	} // if numprocs>1

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "rankid = " << rankid << ", rankA11 = " << rankA11
		<< ", rankA00 = " << rankA00 << ", rankA10 = " << rankA10
		<< ", rankA01 = " << rankA01 << endl;
	}
#endif

	double a11 = 0;
	double a00 = 0;
	double a10 = 0;
	double a01 = 0;
	MPI_Allreduce(&rankA11, &a11, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA00, &a00, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA10, &a10, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&rankA01, &a01, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double randIndex = (a11 + a00) / (a11 + a10 + a01 + a00);
	double m = a11 + a10 + a01 + a00;
	double adjustRandIndex = (a11 - ((a11 + a10) * (a11 + a01)) / m)
			/ ((a11 + a10 + a11 + a01) / 2 - ((a11 + a10) * (a11 + a01)) / m);
	double jaccardIndex = a11 / (a11 + a10 + a01);

#ifdef MPI_DEBUG
	if (OUTPUT_RANK == rankid) {
		cout << "RI = " << randIndex << ", ARI = " << adjustRandIndex
		<< ", JI = " << jaccardIndex << endl;
	}
#endif

	metrics[0] = randIndex;
	metrics[1] = adjustRandIndex;
	metrics[2] = jaccardIndex;
}
//--------------------------------------------------------------------------------------//

/**
 * Get the number of intersected elements between two sets
 */
template<class T>
int MPIMetric<T>::getIntersectionNum(unordered_set<T>& scomm,
		unordered_set<T>& gcomm) {
	int num = 0;
	if (scomm.empty() || gcomm.empty()) {
		return num;
	}

	int scommSize = scomm.size();
	int gcommSize = gcomm.size();

// Miles: Change to HashSet to save time
	unordered_set<T> tmp1;
	unordered_set<T> tmp2;
// Always iterate the set with smaller size
	if (scommSize < gcommSize) {
		tmp1 = scomm;
		tmp2 = gcomm;
	} else {
		tmp1 = gcomm;
		tmp2 = scomm;
	}

	typename unordered_set<T>::const_iterator iter;
	for (iter = tmp1.begin(); iter != tmp1.end(); ++iter) {
		int nodeId = *iter;
		typename unordered_set<T>::const_iterator findIter = tmp2.find(nodeId);
		if (findIter != tmp2.end()) {
			++num;
		}
	}

	return num;
}

/* To allocated a 2d matrix that is continuous in memory */
template<class T>
T** MPIMetric<T>::malloc2d(int row, int column) {
	int size = sizeof(T);
	T **arr;
	arr = (T **) malloc(sizeof(T *) * row + size * row * column);
	if (arr != NULL) {
		T *head;
		head = (T *) arr + sizeof(T *) * row;
		// memset(arr, 0, sizeof(T *) * row + size * row * column);
		while (row--) {
			arr[row] = head + size * row * column;
		}
	}

	return arr;
}

/* Free the memory of 2d matrix */
template<class T>
void MPIMetric<T>::free2d(T **arr) {
	if (arr != NULL) {
		free(arr);
	}
}

template<class T>
MPIMetric<T>::~MPIMetric(void) {

}

#endif /* MPIMETRIC_CPP_ */
