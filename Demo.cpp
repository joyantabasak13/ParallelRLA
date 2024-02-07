#include <iostream>
#include <mpi.h>
#include <sstream>
#include <set>
#include <string>
#include <sys/time.h>
#include <tuple>
#include <time.h>
#include <unistd.h>
#include <utility>
#include <vector>
#include <mutex> 

using namespace std;

int corerank, nprocs, thread_id, numThreads, cxx_procs;



int main(int args, char *argv[]) {
    MPI_Init(&args, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &corerank);
	
	cout<< "Hello from " << corerank << endl; 

    MPI_Finalize();
    return 0;
}