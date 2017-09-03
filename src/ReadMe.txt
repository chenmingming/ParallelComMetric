Authors: Mingming Chen (email: mileschen2008@gmail.com) and Sisi Liu (email: liusisiapply@gmail.com)
Collaborator: Boleslaw K. Szymanski (email: szymab@rpi.edu).

Please cite our paper below if publishing a paper using this toolkit.
Mingming Chen, Sisi Liu, and Boleslaw Szymanski, “Parallel Toolkit for Measuring the Quality of Network Community Structure”, The First European Network Intelligence Conference (ENIC), Wroclaw, Poland, September, 2014, pp. 22-29.



-------------------------------------1. Compile-----------------------------------------------
To compile parallel MPI programs:
Linux: mpic++ -g MPIMetricMain.cpp -o mpimetric
Blue Gene/Q: mpixlcxx -qlanglvl=variadictemplates -o3 MPIMetricMain.cpp -o  mpimetric

To compile parallel Pthreads programs:
Linux: g++ -pthread -g PthreadMetricMain.cpp -o pthreadmetric



---------------------------------------2. Run-------------------------------------------------
To run parallel MPI programs for calculating the metrics with ground truth community structure:
Linux: mpirun -np 4 ./mpimetric metricType realCommunityFile detectedCommunityFile
Blue Gene/Q: srun ./mpimetric metricType realCommunityFile detectedCommunityFile

To run parallel MPI programs for calculating the metrics without ground truth community structure:
Linux: mpirun -np 4 ./mpimetric metricType detectedCommunityFile networkFile [isUnweighted] [isUndirected]
Blue Gene/Q: srun ./mpimetric metricType detectedCommunityFile networkFile [isUnweighted] [isUndirected]


To run parallel Pthreads programs for calculating the metrics with ground truth community structure:
Linux: ./pthreadmetric numThreads metricType realCommunityFile detectedCommunityFile

To run parallel Pthreads programs for calculating the metrics without ground truth community structure:
Linux: ./pthreadmetric numThreads metricType detectedCommunityFile networkFile [isUnweighted] [isUndirected]


Parameters introduction: 
@metricType: metricType=1 for metrics with ground truth community structure; metricType=0 for metrics without ground truth community structure.
@realCommunityFile: it is the file of the ground truth community structure.
@detectedCommunityFile: it is the file of the discovered community structure.
@networkFile: the file of the network.
@isUnweighted: it is optional and default value is 1; isUnweighted=1 for unweighted nework; isUnweighted=0 for weighted network.
@isUndirected: it is optional and default value is 1; isUndirected=1 for undirected network; isUndirected=0 for directed network.
@numThreads: the number of threads adopted in the parallel Pthreads program; Its value should be equal to or larger than 1.



---------------------------------------3. Examples-----------------------------------------------
Example 1: 
mpirun -np 4 ./mpimetric 1 ./dataset/football_true_community.groups ./dataset/football_detected_community.groups

Example 2: 
mpirun -np 4 ./mpimetric 0 ./dataset/football_detected_community.groups ./dataset/football_network.pairs 1 1

Example 3: 
./pthreadmetric 4 1 ./dataset/football_true_community.groups ./dataset/football_detected_community.groups

Example 4:
./pthreadmetric 4 0 ./dataset/football_detected_community.groups ./dataset/football_network.pairs 1 1


-------------------------------------4. Input format----------------------------------------------
(1) The network should be a list of tab or space delimited edges/links with format like "srcNodeId dstNodeId weight" or "srcNodeId\tdstNodeId\tweight".
(2) The ground truth community or detected community file have the format that each line in the file represents a community with node separated with a single space.


-------------------------------------5. Output format---------------------------------------------
(1) Output format of parallel MPI programs for calculating the metrics with ground truth community structure
numProcs total_running_time_information_theory_metrics computation_time msg_passing_time total_running_time_cluster_matching_metrics computation_time msg_passing_time total_running_time_pair_counting_metrics computation_time msg_passing_time
numProcs VI NMI F-measure NVD RI ARI JI

(2) Output format of parallel MPI programs for calculating the metrics without ground truth community structure
numProcs total_running_time 
numProcs modularity modularity_density #intra-edges intra-density contraction #inter-edges expansion conductance

(3) Output format of parallel Pthreads programs for calculating the metrics with ground truth community structure
numThreads total_running_time_information_theory_metrics total_running_time_cluster_matching_metrics total_running_time_pair_counting_metrics
numThreads VI NMI F-measure NVD RI ARI JI

(4) Output format of parallel Pthreads programs for calculating the metrics without ground truth community structure
numThreads total_running_time 
numThreads modularity modularity_density #intra-edges intra-density contraction #inter-edges expansion conductance


-----------------------------------------6. Note--------------------------------------------------
Note: when running the programs, remember to change the predefined clock rate variables in MPIMetricMain.cpp and PthreadMetricMain.cpp to be the value of the clock rate of your own machine. There are used to calculate the running time of the parallel programs.
On shared memory machines, recommend to use parallel Pthreads programs.
