mpirun -np 1 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
mpirun -np 2 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
mpirun -np 4 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
mpirun -np 8 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
mpirun -np 16 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
mpirun -np 24 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
#mpirun -np 64 ./mpimetric 1 ./dataset/LFR/LFR_n100000.txt ./dataset/LFR/SLPA_LFR_n100000.icpm
