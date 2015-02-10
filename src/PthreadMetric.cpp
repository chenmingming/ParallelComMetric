/*
 * PthreadMetric.cpp
 *
 *  Created on: Apr 26, 2014
 *      Author: Sisi
 */

#ifndef PTHREADMETRIC_CPP_
#define PTHREADMETRIC_CPP_

#include "PthreadMetric.h"

template<class T>
PthreadMetric<T>::PthreadMetric() {
	this->numThreads = 1;
	this->isUnweighted = true;
	this->isUndirected = true;

	this->networkFile = "";
	this->realCommunityFile = "";
	this->disCommunityFile = "";

	this->numNodes = 0;
	this->totalWeight = 0;
	this->modularity = 0;
	this->Qds = 0;
	this->intraEdges = 0;
	this->intraDensity = 0;
	this->contraction = 0;
	this->interEdges = 0;
	this->expansion = 0;
	this->conductance = 0;
	this->xentropy = 0;
	this->yentropy = 0;
	this->VI = 0;
	this->NMI = 0;
	this->fMeasure = 0;
	this->NVD = 0;
	this->a11 = 0;
	this->a00 = 0;
	this->a10 = 0;
	this->a01 = 0;
	this->RI = 0;
	this->ARI = 0;
	this->JI = 0;
}

template<class T>
PthreadMetric<T>::PthreadMetric(int numThreads, string networkFile,
		bool isUnweighted, bool isUndirected, string disCommunityFile) {
	this->numThreads = numThreads;
	this->isUnweighted = isUnweighted;
	this->isUndirected = isUndirected;

	this->numNodes = 0;
	this->totalWeight = 0;
	this->modularity = 0;
	this->Qds = 0;
	this->intraEdges = 0;
	this->intraDensity = 0;
	this->contraction = 0;
	this->interEdges = 0;
	this->expansion = 0;
	this->conductance = 0;
	this->xentropy = 0;
	this->yentropy = 0;
	this->VI = 0;
	this->NMI = 0;
	this->fMeasure = 0;
	this->NVD = 0;
	this->a11 = 0;
	this->a00 = 0;
	this->a10 = 0;
	this->a01 = 0;
	this->RI = 0;
	this->ARI = 0;
	this->JI = 0;

	this->networkFile = networkFile;
	this->realCommunityFile = "";
	this->disCommunityFile = disCommunityFile;
}

template<class T>
PthreadMetric<T>::PthreadMetric(int numThreads, string realCommunityFile,
		string disCommunityFile) {
	this->numThreads = numThreads;
	this->isUnweighted = true;
	this->isUndirected = true;

	this->numNodes = 0;
	this->totalWeight = 0;
	this->modularity = 0;
	this->Qds = 0;
	this->intraEdges = 0;
	this->intraDensity = 0;
	this->contraction = 0;
	this->interEdges = 0;
	this->expansion = 0;
	this->conductance = 0;
	this->xentropy = 0;
	this->yentropy = 0;
	this->VI = 0;
	this->NMI = 0;
	this->fMeasure = 0;
	this->NVD = 0;
	this->a11 = 0;
	this->a00 = 0;
	this->a10 = 0;
	this->a01 = 0;
	this->RI = 0;
	this->ARI = 0;
	this->JI = 0;

	this->networkFile = "";
	this->realCommunityFile = realCommunityFile;
	this->disCommunityFile = disCommunityFile;
}

template<class T>
double PthreadMetric<T>::getModularity() {
	return this->modularity;
}

template<class T>
double PthreadMetric<T>::getQds() {
	return this->Qds;
}

template<class T>
double PthreadMetric<T>::getIntraEdges() {
	return this->intraEdges;
}

template<class T>
double PthreadMetric<T>::getIntraDensity() {
	return this->intraDensity;
}

template<class T>
double PthreadMetric<T>::getContraction() {
	return this->contraction;
}

template<class T>
double PthreadMetric<T>::getInterEdges() {
	return this->interEdges;
}

template<class T>
double PthreadMetric<T>::getExpansion() {
	return this->expansion;
}

template<class T>
double PthreadMetric<T>::getConductance() {
	return this->conductance;
}

template<class T>
double PthreadMetric<T>::getVI() {
	return this->VI;
}

template<class T>
double PthreadMetric<T>::getNMI() {
	return this->NMI;
}

template<class T>
double PthreadMetric<T>::getFMeasure() {
	return this->fMeasure;
}

template<class T>
double PthreadMetric<T>::getNVD() {
	return this->NVD;
}

template<class T>
double PthreadMetric<T>::getRI() {
	return this->RI;
}

template<class T>
double PthreadMetric<T>::getARI() {
	return this->ARI;
}

template<class T>
double PthreadMetric<T>::getJI() {
	return this->JI;
}

template<class T>
void* PthreadMetric<T>::metricWithoutGroundTruthAdaptor(void *arg) {
	struct thread_arg<T> *tmp_pointer = (struct thread_arg<T> *) arg;
	return tmp_pointer->classPointer->subWithoutGroundTruthMetricCalculation(
			tmp_pointer);
}

template<class T>
void* PthreadMetric<T>::subWithoutGroundTruthMetricCalculation(
		struct thread_arg<T> *arg) {
	unsigned int threadId = arg->threadId;
	int disComNum = this->disVecCommunities.size();
	// unordered_map<T, int> nodeCommunities;
	unordered_map<int, double> comWeights;
	unordered_map<int, double> comEdges;

	double threadModularity = 0;
	double threadQds = 0;
	double threadIntraEdges = 0;
	double threadIntraDensity = 0;
	double threadContraction = 0;
	double threadInterEdges = 0;
	//double interDensity = 0;
	double threadExpansion = 0;
	double threadConductance = 0;

	// Calculate metrics for this thread
	for (int i = 0; i < disComNum; ++i) {
		if (i % this->numThreads == threadId) {
			comWeights.clear();
			comEdges.clear();
			unordered_set<T> disCommunity = disVecCommunities[i];
			int comId = i;
			double comSize = disCommunity.size();
			double out_incoming_weights = 0;
			typename unordered_set<T>::const_iterator comNodeIt;
			for (comNodeIt = disCommunity.begin();
					comNodeIt != disCommunity.end(); ++comNodeIt) {
				int nodeId = *comNodeIt;
				unordered_map<T, double> nbs = this->network[nodeId];
				typename unordered_map<T, double>::const_iterator nbIt;
				for (nbIt = nbs.begin(); nbIt != nbs.end(); ++nbIt) {
					int nbNodeId = nbIt->first;
					int nbComId = this->disMapCommunities[nbNodeId];
					double weight = nbIt->second;
					double comWeight = comWeights[nbComId];
					comWeights[nbComId] = comWeight + weight;
					double comEdge = comEdges[nbComId];
					comEdges[nbComId] = comEdge + 1;
				}

				if (!this->isUndirected) {
					nbs = this->inNet[nodeId];
					for (nbIt = nbs.begin(); nbIt != nbs.end(); ++nbIt) {
						int nbNodeId = nbIt->first;
						int nbComId = this->disMapCommunities[nbNodeId];

						if (comId != nbComId) {
							double weight = nbIt->second;
							out_incoming_weights += weight;
						}
					}
				}
			} // community node for

			double inWeights = comWeights[comId];
			double inEdges = comEdges[comId];
			double outWeights = 0;
			double splitPenalty = 0;
			typename unordered_map<int, double>::const_iterator weightIter;
			for (weightIter = comWeights.begin();
					weightIter != comWeights.end(); ++weightIter) {
				int nbComId = weightIter->first;
				int nbComSize = this->disVecCommunities[nbComId].size();

				if (comId != nbComId) {
					double comWeight = weightIter->second;
					outWeights += comWeight;
					double comEdge = comEdges[nbComId];
					double sp = (comWeight / this->totalWeight)
							* (comEdge / (comSize * nbComSize));
					splitPenalty += sp;
				}
			}

			// modularity
			if (this->isUndirected) {
				threadModularity += inWeights / this->totalWeight
						- pow((inWeights + outWeights) / this->totalWeight, 2);
			} else {
				threadModularity += inWeights / this->totalWeight
						- ((inWeights + outWeights)
								* (inWeights + out_incoming_weights))
								/ pow(this->totalWeight, 2);
			}

			// Modularity Density Qds
			double inDensity = 0;
			if (comSize > 1) {
				inDensity = inEdges / (comSize * (comSize - 1));
			}
			if (this->isUndirected) {
				threadQds += (inWeights / this->totalWeight) * inDensity
						- pow(
								((inWeights + outWeights) / this->totalWeight)
										* inDensity, 2) - splitPenalty;
			} else {
				threadQds += (inWeights / this->totalWeight) * inDensity
						- (((inWeights + outWeights)
								* (inWeights + out_incoming_weights))
								/ pow(this->totalWeight, 2)) * pow(inDensity, 2)
						- splitPenalty;
			}

			// intra-edges
			if (isUndirected) {
				threadIntraEdges += inWeights / 2;
			} else {
				threadIntraEdges += inWeights;
			}
			// contraction: average degree
			if (inWeights == 0 || comSize == 0) {
				threadContraction += 0;
			} else {
				threadContraction += inWeights / comSize;
			}
			// intra-density
			threadIntraDensity += inDensity;

			threadInterEdges += outWeights;
			// inter-density
			// if (numNodes == srcComSize) {
			// interDensity += 0;
			// } else {
			// interDensity += outWeights
			// / (srcComSize * (numNodes - srcComSize));
			// }
			if (outWeights == 0 || comSize == 0) {
				threadExpansion += 0;
			} else {
				threadExpansion += outWeights / comSize;
			}

			// Avoid that totalInterEdges==0 and communityEdges[i][i]==0
			if (outWeights == 0) {
				threadConductance += 0;
			} else {
				threadConductance += outWeights / (inWeights + outWeights);
			}
		} // end if
	} // end community for

	// Do not update the global variables in threads, must synchronization if doing so
//	this->modularity += threadModularity;
//	this->Qds += threadQds;
//	this->intraEdges += threadIntraEdges;
//	this->intraDensity += threadIntraDensity;
//	this->contraction += threadContraction;
//	this->interEdges += threadInterEdges;
//	this->expansion += threadExpansion;
//	this->conductance += threadConductance;
//	pthread_exit((void*) 0);

	double* results = new double[8];
	results[0] = threadModularity;
	results[1] = threadQds;
	results[2] = threadIntraEdges;
	results[3] = threadIntraDensity;
	results[4] = threadContraction;
	results[5] = threadInterEdges;
	results[6] = threadExpansion;
	results[7] = threadConductance;
	return (void *) results;
}

template<class T>
double PthreadMetric<T>::computeMetricWithoutGroundTruth() {
	// First, clear it in case that it has content.
	this->network.clear();
	this->inNet.clear();
	this->disVecCommunities.clear();
	this->disMapCommunities.clear();
	Reader<T> reader(this->networkFile);
	this->totalWeight = reader.getNetwork(this->isUnweighted,
			this->isUndirected, this->network);
	if (!this->isUndirected) {
		reader.getReversedNetwork(this->isUnweighted, this->isUndirected,
				this->inNet);
	}
	reader.setFileName(this->disCommunityFile);
	this->numNodes = reader.getCommunity(this->disVecCommunities);
	reader.getMapCommunity(this->disMapCommunities);

	double startExeCycle = rdtsc();
	/* pthread_t subThreads[numThreads]; */
	pthread_t *subThreads = new pthread_t[this->numThreads];
	struct thread_arg<T> *thread_args =
			new struct thread_arg<T> [this->numThreads];
	/*
	 Create thread attribute to specify that the main thread needs
	 to join with the threads it creates.
	 */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void *thread_status = NULL;
	int returnCode;

	for (int i = 0; i < this->numThreads; ++i) {
		thread_args[i].threadId = i;
		thread_args[i].classPointer = this;
		returnCode = pthread_create(&subThreads[i], &attr,
				metricWithoutGroundTruthAdaptor, (void *) &thread_args[i]);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr
					<< "computeMetricWithoutGroundTruth::pthread_create return code is "
					<< returnCode << endl;
		}
	}

	/* Wait until the threads finish their task */
	for (int i = 0; i < this->numThreads; ++i) {
		returnCode = pthread_join(subThreads[i], &thread_status);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr
					<< "computeMetricWithoutGroundTruth::pthread_join return code is "
					<< returnCode << endl;
		} else {
			// if success, update the metrics values
			double* results = (double *) thread_status;
			this->modularity += results[0];
			this->Qds += results[1];
			this->intraEdges += results[2];
			this->intraDensity += results[3];
			this->contraction += results[4];
			this->interEdges += results[5];
			this->expansion += results[6];
			this->conductance += results[7];
		}

		free(thread_status);
	}

	pthread_attr_destroy(&attr);
	free(subThreads);
	free(thread_args);

	// Get the final values of metrics
	int numComs = this->disVecCommunities.size();
	this->intraEdges = this->intraEdges / numComs;
	this->intraDensity = this->intraDensity / numComs;
	this->contraction = this->contraction / numComs;
	this->interEdges = this->interEdges / numComs;
	this->expansion = this->expansion / numComs;
	this->conductance = this->conductance / numComs;

	double endExeCycle = rdtsc();
	double exeCycle = endExeCycle - startExeCycle;
	return exeCycle;
}

template<class T>
void* PthreadMetric<T>::infoEntropyMetricAdaptor(void *arg) {
	struct thread_arg<T> *tmp_pointer = (struct thread_arg<T> *) arg;
	return tmp_pointer->classPointer->subInfoEntropyMetricCalculation(
			tmp_pointer);

	// return NULL;
}

template<class T>
void* PthreadMetric<T>::subInfoEntropyMetricCalculation(
		struct thread_arg<T> *arg) {
	//int *tmp_pointer = (int *) arg;
	//int threadId = *tmp_pointer;
	unsigned int threadId = arg->threadId;
	int realComNum = this->realVecCommunities.size();
	int disComNum = this->disVecCommunities.size();
	double threadXentropy = 0;
	double threadYentropy = 0;
	double threadVI = 0;
	double threadNMI = 0;

	// Calculate VI and NMI for this thread
	for (int i = 0; i < realComNum; ++i) {
		if (i % this->numThreads == threadId) {
			unordered_set<T> realCommunity = realVecCommunities[i];
			double ni = realCommunity.size();

			// xentropy
			double tmpEntropy = -(ni / this->numNodes)
					* log2(ni / this->numNodes);
			threadXentropy += tmpEntropy;

			for (int j = 0; j < disComNum; ++j) {
				unordered_set<T> disCommunity = disVecCommunities[j];
				double nj = disCommunity.size();
				double nij = this->getIntersectionNum(realCommunity,
						disCommunity);

				// VI and NMI
				if (nij != 0) {
					double tmpVI = nij * log2((nij * nij) / (ni * nj));
					threadVI += tmpVI;
					double tmpNMI = (nij / this->numNodes)
							* log2((nij * this->numNodes) / (ni * nj));
					threadNMI += tmpNMI;
				}
			}
		}
	}

	// yentropy
	for (int j = 0; j < disComNum; ++j) {
		if (j % this->numThreads == threadId) {
			unordered_set<T> disCommunity = disVecCommunities[j];
			double nj = disCommunity.size();
			double tmpEntropy = -(nj / this->numNodes)
					* log2(nj / this->numNodes);
			threadYentropy += tmpEntropy;
		}
	}

	// Do not update the global variables in threads, must synchronization if doing so
//	this->xentropy += threadXentropy;
//	this->yentropy += threadYentropy;
//	this->VI += threadVI;
//	this->NMI += threadNMI;
//	pthread_exit((void*) 0);

	double* results = new double[4];
	results[0] = threadXentropy;
	results[1] = threadYentropy;
	results[2] = threadVI;
	results[3] = threadNMI;

	return (void *) results;
}

template<class T>
double PthreadMetric<T>::computeInfoEntropyMetric() {
	// First, clear it in case that it has content.
	this->realVecCommunities.clear();
	this->disVecCommunities.clear();
	Reader<T> reader(this->realCommunityFile);
	this->numNodes = reader.getCommunity(realVecCommunities);
	reader.setFileName(this->disCommunityFile);
	reader.getCommunity(disVecCommunities);

	double startExeCycle = rdtsc();
	/* pthread_t subThreads[numThreads]; */
	pthread_t *subThreads = new pthread_t[this->numThreads];
	struct thread_arg<T> *thread_args =
			new struct thread_arg<T> [this->numThreads];
	/*
	 Create thread attribute to specify that the main thread needs
	 to join with the threads it creates.
	 */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void *thread_status = NULL;
	int returnCode;

	for (int i = 0; i < this->numThreads; ++i) {
		thread_args[i].threadId = i;
		thread_args[i].classPointer = this;
		returnCode = pthread_create(&subThreads[i], &attr,
				infoEntropyMetricAdaptor, (void *) &thread_args[i]);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr << "computeInfoEntropyMetric::pthread_create return code is "
					<< returnCode << endl;
		}
	}

	/* Wait until the threads finish their task */
	for (int i = 0; i < this->numThreads; ++i) {
		returnCode = pthread_join(subThreads[i], &thread_status);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr << "computeInfoEntropyMetric::pthread_join return code is "
					<< returnCode << endl;
		} else {
			double* results = (double *) thread_status;
			this->xentropy += results[0];
			this->yentropy += results[1];
			this->VI += results[2];
			this->NMI += results[3];
		}

		free(thread_status);
	}

	pthread_attr_destroy(&attr);
	free(subThreads);
	free(thread_args);

	// Get the final value of VI and NMI
	this->VI = -this->VI / numNodes;
	this->NMI = (2 * this->NMI) / (this->xentropy + this->yentropy);

	double endExeCycle = rdtsc();
	double exeCycle = endExeCycle - startExeCycle;
	return exeCycle;
}

template<class T>
void* PthreadMetric<T>::clusterMatchingMetricAdaptor(void *arg) {
	struct thread_arg<T> *tmp_pointer = (struct thread_arg<T> *) arg;
	return tmp_pointer->classPointer->subClusterMatchingMetricCalculation(
			tmp_pointer);

	// return NULL;
}

template<class T>
void* PthreadMetric<T>::subClusterMatchingMetricCalculation(
		struct thread_arg<T> *arg) {
	unsigned int threadId = arg->threadId;
	int realComNum = this->realVecCommunities.size();
	int disComNum = this->disVecCommunities.size();
	double threadFMeasure = 0;
	double threadNVD = 0;

	// Calculate the first part of fMeasure and NVD of this thread
	for (int i = 0; i < realComNum; ++i) {
		if (i % this->numThreads == threadId) {
			unordered_set<T> realCommunity = this->realVecCommunities[i];
			double ni = realCommunity.size();
			double maxCommon = -1;
			double maxFMeasure = -1;

			for (int j = 0; j < disComNum; ++j) {
				unordered_set<T> disCommunity = this->disVecCommunities[j];
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

			threadNVD += maxCommon;
			double tmpFMeasure = ni * maxFMeasure;
			threadFMeasure += tmpFMeasure;
		} // end if
	} // end for

	// Calculate the second part of fMeasure and NVD of this thread
	for (int i = 0; i < disComNum; ++i) {
		if (i % this->numThreads == threadId) {
			unordered_set<T> disCommunity = this->disVecCommunities[i];
			double maxCommon = -1;

			for (int j = 0; j < realComNum; ++j) {
				unordered_set<T> realCommunity = this->realVecCommunities[j];
				double nij = this->getIntersectionNum(disCommunity,
						realCommunity);

				if (nij > maxCommon) {
					maxCommon = nij;
				}
			}

			threadNVD += maxCommon;
		} // end if
	} // end for

	// Do not update the global variables in threads, must synchronization if doing so
//	this->fMeasure += threadFMeasure;
//	this->NVD += threadNVD;
//	pthread_exit((void*) 0);

	double* results = new double[2];
	results[0] = threadFMeasure;
	results[1] = threadNVD;

	return (void *) results;
}

template<class T>
double PthreadMetric<T>::computeClusterMatchingMetric() {
	// First, clear it in case that it has content.
	this->realVecCommunities.clear();
	this->disVecCommunities.clear();
	Reader<T> reader(this->realCommunityFile);
	this->numNodes = reader.getCommunity(realVecCommunities);
	reader.setFileName(this->disCommunityFile);
	reader.getCommunity(disVecCommunities);

	double startExeCycle = rdtsc();
	/* pthread_t subThreads[numThreads]; */
	pthread_t *subThreads = new pthread_t[this->numThreads];
	struct thread_arg<T> *thread_args =
			new struct thread_arg<T> [this->numThreads];
	/*
	 Create thread attribute to specify that the main thread needs
	 to join with the threads it creates.
	 */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void *thread_status = NULL;
	int returnCode;

	for (int i = 0; i < this->numThreads; ++i) {
		thread_args[i].threadId = i;
		thread_args[i].classPointer = this;
		returnCode = pthread_create(&subThreads[i], &attr,
				clusterMatchingMetricAdaptor, (void *) &thread_args[i]);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr
					<< "computeClusterMatchingMetric::pthread_create return code is "
					<< returnCode << endl;
		}
	}

	/* Wait until the threads finish their task */
	for (int i = 0; i < this->numThreads; ++i) {
		returnCode = pthread_join(subThreads[i], &thread_status);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr << "computeClusterMatchingMetric::pthread_join return code is "
					<< returnCode << endl;
		} else {
			double* results = (double *) thread_status;
			this->fMeasure += results[0];
			this->NVD += results[1];
		}

		free(thread_status);
	}

	pthread_attr_destroy(&attr);
	free(subThreads);
	free(thread_args);

	// Get the final value of fMeasure and NVD
	this->fMeasure = this->fMeasure / this->numNodes;
	this->NVD = 1 - this->NVD / (2 * this->numNodes);

	double endExeCycle = rdtsc();
	double exeCycle = endExeCycle - startExeCycle;
	return exeCycle;
}

template<class T>
void* PthreadMetric<T>::indexMetricAdaptor(void *arg) {
	struct thread_arg<T> *tmp_pointer = (struct thread_arg<T> *) arg;
	return tmp_pointer->classPointer->subIndexMetricCalculation(tmp_pointer);

	// return NULL;
}

template<class T>
void* PthreadMetric<T>::subIndexMetricCalculation(struct thread_arg<T> *arg) {
	unsigned int threadId = arg->threadId;
	vector<T> nodes;
	typename unordered_map<T, int>::const_iterator nodeIt;
	for (nodeIt = this->realMapCommunities.begin();
			nodeIt != this->realMapCommunities.end(); ++nodeIt) {
		nodes.push_back(nodeIt->first);
	}

	double threadA11 = 0;
	double threadA00 = 0;
	double threadA10 = 0;
	double threadA01 = 0;

	// Compute a11, a00, a10, and a01 of this thread
	for (int i = 0; i < this->numNodes; ++i) {
		if (i % this->numThreads == threadId) {
			int inode = nodes[i];
			for (int j = 0; j < this->numNodes; ++j) {
				int jnode = nodes[j];
				if (inode != jnode) {
					int icom = this->disMapCommunities[inode];
					int jcom = this->disMapCommunities[jnode];
					bool dflag = true;
					if (icom != jcom) {
						dflag = false;
					}
					icom = this->realMapCommunities[inode];
					jcom = this->realMapCommunities[jnode];
					bool rflag = true;
					if (icom != jcom) {
						rflag = false;
					}

					if (dflag && rflag) {
						threadA11 += 1;
					} else if (rflag && !dflag) {
						threadA10 += 1;
					} else if (!rflag && dflag) {
						threadA01 += 1;
					} else if (!rflag && !dflag) {
						threadA00 += 1;
					}
				}
			} // inner for
		} // end if
	} // outer for

//	this->a11 += threadA11;
//	this->a00 += threadA00;
//	this->a10 += threadA10;
//	this->a01 += threadA01;
//	pthread_exit((void*) 0);

	double* results = new double[4];
	results[0] = threadA11;
	results[1] = threadA00;
	results[2] = threadA10;
	results[3] = threadA01;

	return (void *) results;
}

template<class T>
double PthreadMetric<T>::computeIndexMetric() {
	// First, clear it in case that it has content.
	this->realMapCommunities.clear();
	this->disMapCommunities.clear();
	Reader<T> reader(this->realCommunityFile);
	this->numNodes = reader.getMapCommunity(this->realMapCommunities);
	reader.setFileName(this->disCommunityFile);
	reader.getMapCommunity(this->disMapCommunities);

	double startExeCycle = rdtsc();
	/* pthread_t subThreads[numThreads]; */
	pthread_t *subThreads = new pthread_t[this->numThreads];
	struct thread_arg<T> *thread_args =
			new struct thread_arg<T> [this->numThreads];
	/*
	 Create thread attribute to specify that the main thread needs
	 to join with the threads it creates.
	 */
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	void *thread_status = NULL;
	int returnCode;

	for (int i = 0; i < this->numThreads; ++i) {
		thread_args[i].threadId = i;
		thread_args[i].classPointer = this;
		returnCode = pthread_create(&subThreads[i], &attr, indexMetricAdaptor,
				(void *) &thread_args[i]);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr << "computeIndexMetric::pthread_create return code is "
					<< returnCode << endl;
		}
	}

	/* Wait until the threads finish their task */
	for (int i = 0; i < this->numThreads; ++i) {
		returnCode = pthread_join(subThreads[i], &thread_status);

		// 0 means success, other values mean failure
		if (returnCode) {
			cerr << "computeIndexMetric::pthread_join return code is "
					<< returnCode << endl;
		} else {
			double* results = (double *) thread_status;
			this->a11 += results[0];
			this->a00 += results[1];
			this->a10 += results[2];
			this->a01 += results[3];
		}

		free(thread_status);
	}

	pthread_attr_destroy(&attr);
	free(subThreads);
	free(thread_args);

	// Get RI, ARI, JI
	double m = this->a11 + this->a10 + this->a01 + this->a00;
	this->RI = (this->a11 + this->a00) / m;
	this->ARI = (this->a11
			- ((this->a11 + this->a10) * (this->a11 + this->a01)) / m)
			/ ((this->a11 + this->a10 + this->a11 + this->a01) / 2
					- ((this->a11 + this->a10) * (this->a11 + this->a01)) / m);
	this->JI = this->a11 / (this->a11 + this->a10 + this->a01);

	double endExeCycle = rdtsc();
	double exeCycle = endExeCycle - startExeCycle;
	return exeCycle;
}

/**
 * Get the number of intersected elements between two sets
 */
template<class T>
int PthreadMetric<T>::getIntersectionNum(unordered_set<T>& scomm,
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

template<class T>
PthreadMetric<T>::~PthreadMetric() {

}

#endif /* PTHREADMETRIC_CPP_ */

