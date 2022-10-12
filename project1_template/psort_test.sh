#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1            # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

SCALE=10000
NP=10

cd /home/cscg/dev/CSC4005/project1_template/

./gen ${SCALE} ./test_data/${SCALE}a.in
./ssort ${SCALE} ./test_data/${SCALE}a.in
echo "ssort done"
mpirun -np ${NP} ./psort ${SCALE} ./test_data/${SCALE}a.in
echo "psort done"
./check ${SCALE} ./test_data/${SCALE}a.in.parallel.out
diff ./test_data/${SCALE}a.in.parallel.out ./test_data/${SCALE}a.in.seq.out
