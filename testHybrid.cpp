#include <iostream>
#include <mpi.h>
#include <pthread.h>
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

// main function for threads
void *threadDriver(void* ptr) {
	int threadID = *static_cast<int*>(ptr);
	cout<< "Rank: " << corerank << " Thread: " << threadID << " Running" << endl;
	return 0;
}


int main(int args, char *argv[]) {
    MPI_Init(&args, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &corerank);
	numThreads = 4;
	cout<< " Rank: " << corerank << " nThreads: " << numThreads << endl;
	usleep(1000);
	pthread_t threads[numThreads-1];
    for (int i = 0; i < numThreads-1; i++)
	{
		int threadID = static_cast<int>(i);
		int iret = pthread_create(&threads[threadID], NULL, threadDriver, &threadID);
		usleep(10);
	}

    MPI_Finalize();
    return 0;
}