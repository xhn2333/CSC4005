g++ sequential.cpp -o seq -O2 -std=c++11
mpic++ mpi.cpp -o mpi -std=c++11
g++ pthread.cpp -lpthread -o pthread -O2 -std=c++11
./seq 1000 1000 100
mpirun -np 4 ./mpi 1000 1000 100
./pthread 1000 1000 100 4