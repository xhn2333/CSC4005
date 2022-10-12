#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1            # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

cd /home/cscg/dev/CSC4005/project1_template/
./gen 10000 ./test_data/10000a.in
./ssort 10000 ./test_data/10000a.in
mpirun -np 20 ./psort 10000 ./test_data/10000a.in
./check 10000 ./test_data/10000a.in.parallel.out