#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1            # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)

# SCALE=100000
# NP=4

declare -i i=0
declare -i j=0
declare -i SCALE=100
declare -i NP=4

cd /home/cscg/dev/CSC4005/project1_template/

for ((i=0; i<=4; i++))
do
    SCALE=100
    for ((j=0; j<=3; j++))
    do
        ./gen ${SCALE} ./test_data/${SCALE}a.in
        mpirun -np ${NP} ./psort ${SCALE} ./test_data/${SCALE}a.in
        ./check ${SCALE} ./test_data/${SCALE}a.in.parallel.out
        SCALE=$SCALE*10
    done
    NP=$NP+4
done
