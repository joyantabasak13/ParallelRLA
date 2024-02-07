// Parallel Record Linkage
// Dynamic Load Balancing
// MPI Status Ignore is okay for now, but have to dealt with later

#include "unionFind.h"
#include "recordComparator.h"
#include "blockingMethods.h"
#include "utilities.h"
#include "sorting.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime> 
#include <fstream>
#include <iostream>
#include <map>
#include <mpi.h>
#include <pthread.h>
#include <queue>
#include <set>
#include <stdio.h>
#include <string>
#include <sys/time.h>
#include <random>
#include <time.h>
#include <tuple>
#include <unistd.h>
#include <utility>
#include <vector>


using namespace std;
using namespace boost;

#define ROOT_ID				0
#define UNIONFIND_ARR		1
#define KMER				3
#define WORK_COMPLETE		4
#define WORK_ASSINGED		5
#define SEND_INDICES_TOBE_SORTED 101
#define SEND_STRINGS_TOBE_SORTED 102
#define SEND_SORTED_INDICES 103

// auto rng = std::default_random_engine {};

int alphabetSize = 36;
int attributes = 5;
int blockingAttr = 1;
int compDividerMultiplier = 100;
int coreRank;
int datasizeAsChar;
int lenMax;
int numCores;
int numThreads = 5;
int perCoreStaticAllocatedChunk = 10;
int threshold = 99;
int totalRecords;
int totalUniqueRecords;
int totalUniqueBlocks;
int world;

long long int totalCompRequired;
long long int totalEdges = 0;

char* superString;
int* assignedRecordIndices;
int* uniqueRecordIDs;


map<int, vector<int> > approxConnectedComponents;

vector<pair<int, string> > combinedData;

vector<vector<int> > exactmatches;
vector<vector<int> > finalConnectedComponents;
vector<vector<pair<long long int, int> > > assignedRecordDomain;
vector<vector<pair<long long int, int> > > assignedRecordRange;
vector<vector<string>> uniqueRecords;
vector<vector<string> > vec2D;


UnionFind uf;
RecordComparator recCom;

// Sorting Variables

vector<pair<int, string>> sortedSampledRecords;
vector<pair<int, string>> chosenSortedSampledRecords;
vector<vector<int>> bucketInds;
vector<vector<int>> bucketAssignments;
vector<int> bucketStarts;
vector<int> totalRecsAssigned;
int sampleSize;
vector<vector<pair<int, string> > > sortedStrings;
vector<int> totalCoreRecs;


class ThreadData {
    public:
		int theadID;
		vector<pair<int, string>> sortedStrings_t;
};

class CompThreadData {
    public:
		int theadID;
		UnionFind uf_t;
		BlockingMethods *bm;

};

queue<int> workQ;
pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;	/* mutex lock for buffer */
pthread_cond_t cond = PTHREAD_COND_INITIALIZER; /* consumer waits on this cond var */
pthread_cond_t prod_cond = PTHREAD_COND_INITIALIZER; /* producer waits on this cond var */
static int chunkProcessed =0;

void getUniqueEntriesGlobalMaster() {
	uniqueRecordIDs = new int[totalUniqueRecords];
	uniqueRecords.resize(totalUniqueRecords);
	for(size_t i=0; i<totalUniqueRecords; i++){
		uniqueRecords[i].resize(attributes);
	}

	for (size_t i = 0; i < totalUniqueRecords; i++) {
        uniqueRecordIDs[i] = exactmatches[i][0];
		for(size_t j=0; j<attributes; j++){
			uniqueRecords[i][j] = vec2D[exactmatches[i][0]][j];
		}
    }
}

void getUniqueEntriesWorkers() {
	uniqueRecords.resize(totalUniqueRecords);
	for(size_t i=0; i<totalUniqueRecords; i++){
		uniqueRecords[i].resize(attributes);
	}

	for (size_t i = 0; i < totalUniqueRecords; i++) {
		for(size_t j=0; j<attributes; j++){
			uniqueRecords[i][j] = vec2D[uniqueRecordIDs[i]][j];
		}
    }
}


void getSuperStringSize(){
	datasizeAsChar = 0;
	for (size_t i = 0; i < totalUniqueRecords; i++) {
        for(size_t j = 0; j< attributes; j++) {
            datasizeAsChar += std::strlen(uniqueRecords[i][j].c_str()) +1;
        }
    }
	// cout<< "Unique Records Super String Size: " << datasizeAsChar << endl;
}


void printWorkChunkSizes(){
	long long int totalCompEst = 0;
	long long int curCompEst = 0;
	long long int minChunk = LLONG_MAX;
	long long int maxChunk = 0;
	int totalChunks = (numCores-1)*compDividerMultiplier;
	// cout<< "Total chunks: " << totalChunks << endl;
	for(int c = 0; c < totalChunks; c++){
		curCompEst = 0;
		for(int i = 0; i < assignedRecordDomain[c].size(); i++) {
			int n = assignedRecordDomain[c][i].second;
			curCompEst+= n * assignedRecordRange[c][i].second - ceil(((n-2)*(n-1))/2) ;
		}
		// cout<< "Chunk: " << c << " estimated comp: " << curCompEst << endl;
		totalCompEst += curCompEst;
		minChunk = minChunk>curCompEst?curCompEst:minChunk;
		maxChunk = maxChunk<curCompEst?curCompEst:maxChunk;
	}
	// cout<<"Loop ended" << endl;
	// cout<< "Total Estimated Comp: " << totalCompEst << endl;
	// cout<< "Max Estimated Comp: " << maxChunk << endl;
	// cout<< "Min Estimated Comp: " << minChunk << endl;
	// cout<< "Load Balancing difference: " << maxChunk-minChunk << endl;
}


void doRequestedRecordComp(int requestedChunk, BlockingMethods &bm) {
	for(int i = 0; i < assignedRecordDomain[requestedChunk].size(); i++) {
		for(int j=assignedRecordDomain[requestedChunk][i].first; j < assignedRecordDomain[requestedChunk][i].first + assignedRecordDomain[requestedChunk][i].second; j++) {
			for(int k=j+1; k < assignedRecordRange[requestedChunk][i].first + assignedRecordRange[requestedChunk][i].second; k++){
				int recid_j = bm.blockingIDList[j].second;
				int recid_k = bm.blockingIDList[k].second;
				if(!uf.isConnected(recid_j, recid_k)){
				 	if(recCom.isLinkageOk(uniqueRecords[recid_j], uniqueRecords[recid_k])) {
						uf.weightedUnion(recid_j, recid_k);
				 		totalEdges++;
				 	}
				}
			}
		}
	}
}


void mergeEdges(UnionFind &ufTemp) {
	for(int i=0; i< totalUniqueRecords; i++) {
		int root_i = uf.find(i);
		int root_mt = ufTemp.find(i);
		if( root_i != root_mt){
			uf.weightedUnion(root_i, root_mt);
		}
	}
}


void getBlockingKeyData(vector<pair<int, string> > &blockingKeyData) {
	string strSample(50, '0');
	int max = 0;
	for (int i = 0; i< totalUniqueRecords; i++) {
		pair<int, string> p;
		p.first = i;
		p.second = uniqueRecords[i][blockingAttr];
		blockingKeyData[i] = p;
		if (max<blockingKeyData[i].second.size()) {
			max = blockingKeyData[i].second.size();
		}
	}
	lenMax = max;
	// Padding to make all characters same size
    for (int i = 0; i < totalUniqueRecords; ++i) {
		int lenDiff		= lenMax - blockingKeyData[i].second.length();
		if(lenDiff > 0)
			blockingKeyData[i].second	+= strSample.substr(0, lenDiff);
	}
}


void doSortedComp() {
	clock_t currTS_f0	= clock();
	vector<pair<int,string> > headlessCopies;
	headlessCopies.resize(2*totalUniqueRecords);
	getBlockingKeyData(headlessCopies);
	clock_t currTS_f1	= clock();
	for(int i=0; i< uniqueRecords.size(); i++) {
		headlessCopies[totalUniqueRecords+i].first = i;
		headlessCopies[totalUniqueRecords+i].second = headlessCopies[i].second.substr(1,lenMax) + '0';
	}
	clock_t currTS_f2	= clock();
	radixSort(lenMax, headlessCopies);
	clock_t currTS_f3	= clock();
	int extraEdges = 0;
	for (int i = 1; i < headlessCopies.size(); i++) {
		if (headlessCopies[i-1].second.compare(headlessCopies[i].second) == 0) {
			int recID_i = headlessCopies[i-1].first;
            int recID_j = headlessCopies[i].first;
            if (!uf.isConnected(recID_i, recID_j)) {
                if(recCom.isLinkageOk(uniqueRecords[recID_i], uniqueRecords[recID_j])) {
					uf.weightedUnion(recID_i, recID_j);
					extraEdges++;
					totalEdges++;
				}
			}
		}
	}
	clock_t currTS_f4	= clock();

	double dataRead_t	= (double)(currTS_f1 - currTS_f0) / CLOCKS_PER_SEC;
	double getHeadLessCp= (double)(currTS_f2 - currTS_f1) / CLOCKS_PER_SEC;
	double sorting_t	= (double)(currTS_f3 - currTS_f2) / CLOCKS_PER_SEC;
	double comp_t		= (double)(currTS_f4 - currTS_f3) / CLOCKS_PER_SEC;

	cout<< "DataRead time: " << dataRead_t << endl;
	cout<< "getHeadLessCp time: " << getHeadLessCp << endl;
	cout<< "sorting_t time: " << sorting_t << endl;
	cout<< "comp_t time: " << comp_t << endl;
	cout<< "LenMax: " << lenMax << endl; 
	// cout<< "Edges Added by Sorting: "<< extraEdges << endl;
}


void findConnComp()
{
	int root, edgeTotal;

    for (int i = 0; i < totalUniqueRecords; i++)
	{
		root = uf.find(i);
		approxConnectedComponents[root].push_back(i);
	}

    cout<< " Single Linkage Connected Components: " << approxConnectedComponents.size()<<endl;
}


void findFinalConnectedComp() {
    int totalNodes = 0;
    int pairsAccessed = 0;
    for (auto const& p : approxConnectedComponents) {
        pairsAccessed++;
        int numComponents = p.second.size();
		totalNodes+=numComponents;
        bool distmat[numComponents][numComponents];
        vector<vector<string>> dataArr(numComponents); // to make cache-efficient, keep records in a row

		// Copy Data in cluster
        for(int c=0; c<p.second.size(); c++) {
            dataArr[c] = uniqueRecords[p.second[c]];
        };

		// generate a 2D matrix filled with all pair comparision results
        for (int i =0; i<numComponents; i++) {
            distmat[i][i] = true;
            for (int j = i+1; j < numComponents; j++)
            {
                if (recCom.isLinkageOk(dataArr[i], dataArr[j])) {
                    distmat[i][j] = true;
                    distmat[j][i] = true;
                } else {
                    distmat[i][j] = false;
                    distmat[j][i] = false;
                }
            }
        }

        bool nodesConsidered[numComponents];
        for(int i=0; i<numComponents; i++) {
            nodesConsidered[i] = false;
        }

        for(int i=0; i<numComponents; i++) {
            if(nodesConsidered[i] == false) {
                vector<int> connectedComponent;
                connectedComponent.push_back(p.second[i]);
                nodesConsidered[i] = true;
                for(int j=0; j<numComponents; j++) {
                    if ((distmat[i][j] == true) && (nodesConsidered[j]==false)) {
                        connectedComponent.push_back(p.second[j]);
                        nodesConsidered[j] = true;
                    }
                }
                finalConnectedComponents.push_back(connectedComponent);
            }
        }
    }
    // cout<< "Total Nodes: "<< totalNodes << " unique records: " << totalUniqueRecords << endl;
    cout<< "Complete Linkage Connected Components: "<< finalConnectedComponents.size()<<endl;
}


void writeFinalConnectedComponentToFile(string& result_file_name) {
	ofstream out_file;
    out_file.open(result_file_name);

	for(int i = 0; i< finalConnectedComponents.size(); i++) {
        for(int j = 0; j< finalConnectedComponents[i].size(); j++) {
            for(int k=0; k< exactmatches[finalConnectedComponents[i][j]].size(); k++) {
                out_file << uniqueRecords[finalConnectedComponents[i][j]][0] << ",";
			}
		}
        out_file<< "\n";
	}
	out_file.close();
}


void doNormalThreadedSorting(ThreadData &myTData) {
	int tID = myTData.theadID;
	// vector<pair<int, string> > *sortedStrings_t = myTData.sortedStrings_t;
	// cout<< "CoreRank: "<< coreRank << "Thread: "<< tID << " has " << totalRecsAssigned[coreRank*numThreads + tID] << " records" << endl;
	// Moved memory allocation to allocate thread operations
	myTData.sortedStrings_t.resize(totalRecsAssigned[coreRank*numThreads + tID]);
	int counter = 0;

	for(int i = 0; i< bucketAssignments[coreRank*numThreads + tID].size(); i++) {
		int assignedBucketId = bucketAssignments[coreRank*numThreads + tID][i];
		if((assignedBucketId < 0) || (assignedBucketId>bucketInds.size())){
			cout<< "ERROR" << endl;
		} 
		for(int j=0; j<bucketInds[assignedBucketId].size(); j++){
			myTData.sortedStrings_t[counter] = combinedData[bucketInds[assignedBucketId][j]];
			counter++;
		}
	}
	// cout<< "Counter: " << counter << " size: " << myTData.sortedStrings_t.size() << endl;

	radixSort_std(lenMax, myTData.sortedStrings_t);
}


void *threadDriverSort(void* ptr) {
	ThreadData *myTData = (ThreadData*)(ptr);
	doNormalThreadedSorting(*myTData);
	return 0;
}

void *threadDriverComp(void* ptr) {
	CompThreadData *myTData = (CompThreadData*)(ptr);
	myTData->uf_t.setVariable(totalUniqueRecords);
	while(1){
		pthread_mutex_lock(&mtx);

		while (workQ.empty()){
			int rc = pthread_cond_wait(&cond, &mtx);
		}
		int workC = workQ.front();
		if(workC == -1) {
			pthread_mutex_unlock(&mtx);
			pthread_cond_signal (&cond);
			break;
		}
		workQ.pop();
		cout<< "CoreRank: " << coreRank << " Thread: " << myTData->theadID << " Processed work chunk: " << workC << endl;
		chunkProcessed++;
		pthread_mutex_unlock(&mtx);
		pthread_cond_signal (&cond);
	}
}


void parallelSampleSort_2(){
	// one thread per node won't work but every other thread will work
	// So we must assume there are numThreads + 1 thread working per node in reality
	sampleSize = numThreads*numCores - 1;

	int regionPerSample = 10;
	int effectiveSampleSize = (sampleSize + 1)*regionPerSample;
	int* sampleInd;
	sampleInd = new int[effectiveSampleSize];
	
	if (coreRank == 0) {
		cout<< "SampleSize: " << sampleSize << " Effective Sample Size: " << effectiveSampleSize << endl;
		// Broadcast this samples, others can be calculated deterministically
		generateRandomSample(effectiveSampleSize, sampleInd, combinedData.size());
	}

	MPI_Bcast(sampleInd, effectiveSampleSize, MPI_INT, 0, MPI_COMM_WORLD);

	getRandomSamples(sortedSampledRecords, sampleInd, combinedData, effectiveSampleSize);
	sortRandomSamples(sortedSampledRecords, sampleInd, lenMax);
	for(int i=regionPerSample-1; i< effectiveSampleSize-regionPerSample+1; i++){
		if((i%(regionPerSample-1)) == 0){
			chosenSortedSampledRecords.push_back(sortedSampledRecords[i]);
		}
	}
	// cout<< "Rank: " << coreRank << "ChosenSortedSamples: " << chosenSortedSampledRecords.size() << endl;
	doBinarySearchAndFindBuckets(chosenSortedSampledRecords, combinedData, bucketInds, bucketStarts);
	getBucketAssignment(bucketInds, bucketAssignments, totalRecsAssigned, sampleSize);


	// Sorting will be done by threads
	// Master thread in the node will collect their results
	// Master threads of each node will send the result to the global master thread

	pthread_t threads[numThreads];
	ThreadData tData[numThreads];

	for (int i = 0; i < numThreads; i++)
	{
		int threadID = static_cast<int>(i);
		tData[i].theadID = static_cast<int>(i);
		int iret = pthread_create(&threads[threadID], NULL, threadDriverSort, &tData[i]);
		usleep(10);
	}


	totalCoreRecs.resize(numCores);
	for(int i=0; i<numCores; i++){
		int totRec = 0;
		for(int j=0; j<numThreads; j++){
			totRec += totalRecsAssigned[i*numThreads + j];
		}
		totalCoreRecs[i] = totRec;
	}
	int *partiallySortedIndices;

	if(coreRank != ROOT_ID) {
		partiallySortedIndices = new int[totalCoreRecs[coreRank]];
	} else {
		partiallySortedIndices = new int[combinedData.size()];
	}

	for (int i = 0; i < numThreads; i++)
	{
		int threadID = static_cast<int>(i);
		pthread_join(threads[threadID], NULL);
	}


	int ind = 0;
	for(int i=0; i<numThreads; i++){
		for(int j=0; j<totalRecsAssigned[coreRank*numThreads + i]; j++){
			partiallySortedIndices[ind] = tData[i].sortedStrings_t[j].first;
			ind++;
		}
		// cout<< totalRecsAssigned[coreRank*numThreads + i] << " by " << ind << endl;
	}
	// cout<< "CoreRank: " << coreRank << "Total: " << totalCoreRecs[coreRank] << " total count " << ind << endl; 


	// cout<< "READY TO SEND OUT SORTED STRINGS" << endl;
	if(coreRank != ROOT_ID){
		MPI_Send(partiallySortedIndices, totalCoreRecs[coreRank], MPI_INT, 0, SEND_SORTED_INDICES, MPI_COMM_WORLD);
	}
	if(coreRank == ROOT_ID) {
		for(int i=1; i<numCores; i++){
			int *recIds = new int[totalCoreRecs[i]];
			MPI_Recv(recIds, totalCoreRecs[i], MPI_INT, i, SEND_SORTED_INDICES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			int start = 0;
			for(int j = 0; j<i; j++){
				start += totalCoreRecs[j];
			}
			for(int j=0; j<totalCoreRecs[i]; j++){
				partiallySortedIndices[start+j] = recIds[j];

			}
		}

		vector<pair<int, string>> tempSortedStrings;
		tempSortedStrings.resize(combinedData.size());
		int ind2 =0;
		for(int i = 0; i<bucketAssignments.size(); i++){
			for(int j=0; j<bucketAssignments[i].size(); j++){
				int startInd = bucketStarts[bucketAssignments[i][j]];
				for(int k =0; k<bucketInds[bucketAssignments[i][j]].size(); k++){
					tempSortedStrings[startInd+k] = combinedData[partiallySortedIndices[ind2]];
					ind2++;
				}
			}
		}
		combinedData = tempSortedStrings;
	}
}



int main(int argc, char** argv) {
	// Start timing from here
	clock_t currTS_p0	= clock();
	clock_t currTS_p1;
	clock_t currTS_p2;
	clock_t currTS_p3;
	clock_t currTS_p4;
	clock_t currTS_p5;
	clock_t currTS_p6;

	int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if(provided < MPI_THREAD_FUNNELED)
    {
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    else
    {
        printf("The threading support level corresponds to that demanded.\n");
    }


    MPI_Comm_rank(MPI_COMM_WORLD, &coreRank);
    MPI_Comm_size(MPI_COMM_WORLD, &world);
	// string filePath = "/Users/joyanta/Documents/Research/Record_Linkage/codes/my_codes/ds_single_datasets/firstName_LastName_DS/";
	string filePath = "/home/job22011/RecordLinkage/data/";
	// string filePath = "/gpfs/scratchfs1/sar02010/job22011/RecordLinkage/data/";
	string fileName = argv[1];
	filePath = filePath + argv[1];

	numCores = world;
	vector<int> attrThresholds{1,1,1,1};

	recCom.setComparator(1, attrThresholds);

	getFormattedDataFromCSV(filePath, vec2D, attributes, totalRecords);

	double dataRead_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    cout<< "Rank " << coreRank <<" DataRead time "<< dataRead_t << endl;
	currTS_p1	= clock();
	// concatenate all attributes for a record
	cout<< "Total Records: " << vec2D.size() << endl;
	lenMax = 0;
	totalUniqueRecords = 0;
	getCombinedData(vec2D, combinedData, lenMax);
	double comb_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
	cout<< "Rank " << coreRank <<" Got Combined Data "<< comb_t << endl;

	// Sort the concatenated records
	parallelSampleSort_2();

	if(coreRank == ROOT_ID){
		double total_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
		cout<< "TOTAL TIME: "<< total_t << endl;
	}


    if(coreRank == ROOT_ID){
		double sorting_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
		cout<< "Sorting time "<< sorting_t << endl;
		// Get Unique Records
        getExactMatches(totalUniqueRecords, exactmatches, combinedData);
        getUniqueEntriesGlobalMaster();
		cout<< "Total Unique Records: " << totalUniqueRecords << endl;
        double exactClustering_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
        cout<< "De-duplication Time "<< exactClustering_t << endl;
	}
	currTS_p2	= clock();


	MPI_Bcast(&totalUniqueRecords, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	if(coreRank!= ROOT_ID){
		uniqueRecordIDs = new int[totalUniqueRecords];
	}
	MPI_Bcast(uniqueRecordIDs, totalUniqueRecords, MPI_INT, ROOT_ID, MPI_COMM_WORLD);

	uf.setVariable(totalUniqueRecords);

	if(coreRank!= ROOT_ID){
		getUniqueEntriesWorkers();
	}

	// if(coreRank == ROOT_ID){
	// 	currTS_p3	= clock();
	// 	doSortedComp();
	// 	double sortedComp_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
	// 	cout<< "Rank: " << coreRank << "Superblocking Sorting Time "<< sortedComp_t << endl;
	// }
	
	
	BlockingMethods bm;
	bm.setBlocking(KMER, alphabetSize, blockingAttr, uniqueRecords);
	bm.getSuperBlockingIDArray();
	double blockingArray_t	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
	cout<< "Rank: " << coreRank << " Blocking Time "<< blockingArray_t << endl;

	// Get Sorted Blocking ID Array
	bm.sortBlockingIDArray();
	bm.removeRedundentBlockingID();

	// Find Block assignments
	bm.findBlockBoundaries();
	bm.findBlockAssignments(numCores, numThreads, compDividerMultiplier, assignedRecordDomain, assignedRecordRange);

	// Add threading code here
	pthread_t wthreads[numThreads];
	CompThreadData tData[numThreads];

	for(int i=0; i<perCoreStaticAllocatedChunk*numThreads; i++){
		workQ.push(coreRank*perCoreStaticAllocatedChunk*numThreads + i);
	}

	for (int i = 0; i < numThreads; i++)
	{
		int threadID = static_cast<int>(i);
		tData[i].theadID = static_cast<int>(i);
		tData[i].bm = &bm;
		int iret = pthread_create(&wthreads[threadID], NULL, threadDriverComp, &tData[i]);
		usleep(10);
	}


	// clock_t currTS_p8	= clock();
	// for(int i=0; i<perCoreStaticAllocatedChunk; i++){
	// 	int selectedChunk = (coreRank-1)*perCoreStaticAllocatedChunk + i;
	// 	doRequestedRecordComp(selectedChunk, bm);
	// }

	while(workQ.size()>=numThreads){
		pthread_cond_wait (&cond, &mtx);
	}
	
	if(coreRank!= ROOT_ID){
		int selectedChunk = 0;
		while(selectedChunk != -1){
			MPI_Send(&coreRank, 1, MPI_INT, ROOT_ID, WORK_COMPLETE, MPI_COMM_WORLD);
			MPI_Recv(&selectedChunk, 1, MPI_INT, ROOT_ID, WORK_ASSINGED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			pthread_mutex_lock(&mtx);
			workQ.push(selectedChunk);
			pthread_mutex_unlock(&mtx);
			pthread_cond_signal(&cond, &mtx);
		}
	}

	if(coreRank == ROOT_ID){
		int received=-1;
		int assignedChunk = numCores*numThreads*perCoreStaticAllocatedChunk;
		MPI_Request requests[numCores-1];
		for(int i = 0; i< numCores-1; i++){
			MPI_Irecv(&received, 1, MPI_INT, i+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[i]);
		}
		while(assignedChunk < numCores*numThreads*compDividerMultiplier){
			int index_count;
            int indices[numCores-1];

            MPI_Waitsome(numCores-1, requests, &index_count, indices, MPI_STATUSES_IGNORE);
			pthread_cond_signal(&cond, &mtx);

            for(int i = 0; i < index_count; i++)
            {
				if(assignedChunk < numCores*numThreads*compDividerMultiplier{
					MPI_Send(&assignedChunk, 1, MPI_INT, indices[i]+1, WORK_ASSINGED, MPI_COMM_WORLD);
					MPI_Irecv(&received, 1, MPI_INT, indices[i]+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[indices[i]]);
					assignedChunk++;
				}
            }
			pthread_mutex_lock(&mtx);
			workQ.push(selectedChunk);
			pthread_mutex_unlock(&mtx);
			pthread_cond_signal(&cond, &mtx);

		}
        MPI_Waitall(numCores-1, requests, MPI_STATUSES_IGNORE);
		int endFlag = -1;
		for(int i=1; i<numCores; i++ ){
			MPI_Send(&endFlag, 1, MPI_INT, i, WORK_ASSINGED, MPI_COMM_WORLD);
		}
	}

	for (int i = 0; i < numThreads; i++)
	{
		int threadID = static_cast<int>(i);
		pthread_join(wthreads[threadID], NULL);
	}

	MPI_Finalize();
	return 0;

	if(coreRank != ROOT_ID) {
		MPI_Send(uf.parentInd, totalUniqueRecords, MPI_INT, 0, UNIONFIND_ARR, MPI_COMM_WORLD);
		double comp_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
		cout<< "Rank: " << coreRank << "Comparision Time "<< comp_t << endl;

	}

	if(coreRank == ROOT_ID) {
		currTS_p4	= clock();
		if(world>1) {
			for(int i=1; i<numCores; i++){
				int tempArr[totalUniqueRecords];
				MPI_Recv(tempArr, totalUniqueRecords, MPI_INT, i, UNIONFIND_ARR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				// Merge UnionFind arr
				UnionFind ufTemp;
				ufTemp.setVariable(totalUniqueRecords);
				ufTemp.setParentArr(tempArr);
				mergeEdges(ufTemp);
			}
		}

		findConnComp();
		double singleLinkage_t	= (double)(clock() - currTS_p4) / CLOCKS_PER_SEC;
		cout<< "Single linkage Time "<< singleLinkage_t << endl;

		currTS_p5	= clock();
		findFinalConnectedComp();
	
		double connCompFind_t	= (double)(clock() - currTS_p5) / CLOCKS_PER_SEC;
		cout<< "Complete Linkage Time "<< connCompFind_t << endl;

		currTS_p6	= clock();
		string outname = "Out_" + std::to_string(world) + "_"+ argv[1] + ".csv";
		writeFinalConnectedComponentToFile(outname);
		double fileWriteTime	= (double)(clock() - currTS_p6) / CLOCKS_PER_SEC;
		cout<< "Data Write Time "<< fileWriteTime << endl;
	}

	if(coreRank == ROOT_ID){
		double total_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
		cout<< "TOTAL TIME: "<< total_t << endl;
	}

	MPI_Finalize();
	return 0;
}