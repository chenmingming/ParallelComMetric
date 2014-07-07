#!/bin/sh

sbatch --time=30 --nodes=1 --ntasks=1 --partition=small ./srun_script.sh
sleep 1800

sbatch --time=30 --nodes=1 --ntasks=2 --partition=small ./srun_script.sh
sleep 1800

sbatch --time=30 --nodes=1 --ntasks=4 --partition=small ./srun_script.sh
sleep 1800

sbatch --time=30 --nodes=1 --ntasks=8 --partition=small ./srun_script.sh
sleep 1800

sbatch --time=30 --nodes=1 --ntasks=16 --partition=small ./srun_script.sh
sleep 1200

sbatch --time=30 --nodes=2 --ntasks=32 --partition=small ./srun_script.sh
sleep 800

sbatch --time=30 --nodes=4 --ntasks=64 --partition=small ./srun_script.sh
sleep 500

sbatch --time=30 --nodes=8 --ntasks=128 --partition=small ./srun_script.sh
sleep 400

sbatch --time=30 --nodes=16 --ntasks=256 --partition=small ./srun_script.sh
sleep 200

sbatch --time=30 --nodes=32 --ntasks=512 --partition=small ./srun_script.sh
sleep 100

sbatch --time=30 --nodes=64 --ntasks=1024 --partition=small ./srun_script.sh
sleep 50

sbatch --time=30 --nodes=128 --ntasks=2048 --partition=medium ./srun_script.sh
sleep 50

sbatch --time=30 --nodes=256 --ntasks=4096 --partition=medium ./srun_script.sh
