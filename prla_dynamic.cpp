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

// auto rng = std::default_random_engine {};

int alphabetSize = 36;
int attributes = 5;
int blockingAttr = 1;
int compDividerMultiplier = 100;
int coreRank;
int datasizeAsChar;
int lenMax;
int numCores;
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


void getUniqueEntries() {
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


void getSuperStringSize(){
	datasizeAsChar = 0;
	for (size_t i = 0; i < totalUniqueRecords; i++) {
        for(size_t j = 0; j< attributes; j++) {
            datasizeAsChar += std::strlen(uniqueRecords[i][j].c_str()) +1;
        }
    }
	// cout<< "Unique Records Super String Size: " << datasizeAsChar << endl;
}

// Super string is a single string generated 
// cancatenating all attributes of unique records seperated by a special sign 
void getSuperString(){
	getSuperStringSize();
	superString = new char[datasizeAsChar];
	int ind = 0;
	int totalTokens = 0;
	for (size_t i = 0; i < totalUniqueRecords; i++) {
        for(size_t j = 0; j< attributes; j++) {
			for(size_t k=0; k< std::strlen(uniqueRecords[i][j].c_str()); k++){
				superString[ind++] = uniqueRecords[i][j][k];
			}
			superString[ind++] = '|';
			totalTokens++;
        }
    }
}


void populateUniqueRecordsFromSuperString(){
    string superStringStr(superString);
	uniqueRecords.resize(totalUniqueRecords);
	for(size_t i=0; i<totalUniqueRecords; i++){
		uniqueRecords[i].resize(attributes);
	}
	vector<string> attrStrs;
    boost::split(attrStrs, superStringStr, boost::is_any_of("|"));
	int ind = 0;
	// cout<< "Total tokens: " << attrStrs.size() << endl;
	// split creates an empty token at the end(after last terminating sign)
	for(int i=0; i<attrStrs.size()-1; i++) {
		if((i != 0) && (i%attributes == 0)){
			ind++;
		}
		// cout<< "Index accessing: " << ind << " Token: " << attrStrs[i] << " substring no: " << i <<endl;
		uniqueRecords[ind][i%attributes] = attrStrs[i];
	}
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


int main(int argc, char** argv) {
	// Start timing from here
	clock_t currTS_p0	= clock();
	clock_t currTS_p1;
	clock_t currTS_p2;
	clock_t currTS_p3;
	clock_t currTS_p4;
	clock_t currTS_p5;
	clock_t currTS_p6;


    MPI_Init(NULL, NULL);

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
	// getCombinedData(vec2D, combinedData, lenMax);
	double comb_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
	cout<< "Rank " << coreRank <<" Got Combined Data "<< comb_t << endl;
    // if (coreRank == ROOT_ID) {
	// 	cout<< "Total Cores: " << world << endl;
	// 	// Read data only on root core;
	// 	// Find out unique records and distribute them to other cores
	// 	// getFormattedDataFromCSV(filePath, vec2D, attributes, totalRecords);

    //     double dataRead_t	= (double)(clock() - currTS_p0) / CLOCKS_PER_SEC;
    //     cout<< "DataRead time "<< dataRead_t << endl;

	// 	currTS_p1	= clock();
	// 	// concatenate all attributes for a record
	// 	getCombinedData(vec2D, combinedData, lenMax);
	// }
	// Sort the concatenated records
	// parallelSampleSort_1(lenMax, numCores, combinedData, coreRank);
	// radixSort(lenMax, combinedData);

	MPI_Finalize();
	return 0;

    if(coreRank == ROOT_ID){
		double sorting_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
		cout<< "Sorting time "<< sorting_t << endl;
		// Get Unique Records
        getExactMatches(totalUniqueRecords, exactmatches, combinedData);
        getUniqueEntries();
        double exactClustering_t	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
        cout<< "De-duplication Time "<< exactClustering_t << endl;
		getSuperString();
	}
	currTS_p2	= clock();

	
	MPI_Bcast(&totalUniqueRecords, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&attributes, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	MPI_Bcast(&datasizeAsChar, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
	if(coreRank != ROOT_ID) {
		superString = new char[datasizeAsChar];
	}
	MPI_Bcast(superString, datasizeAsChar, MPI_CHAR, ROOT_ID, MPI_COMM_WORLD);

	uf.setVariable(totalUniqueRecords);

	if(coreRank == ROOT_ID){
		currTS_p3	= clock();
		doSortedComp();
		double sortedComp_t	= (double)(clock() - currTS_p3) / CLOCKS_PER_SEC;
		cout<< "Rank: " << coreRank << "Superblocking Sorting Time "<< sortedComp_t << endl;
	}

	if(coreRank!= ROOT_ID){
		populateUniqueRecordsFromSuperString();
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
		bm.findBlockAssignments(numCores, compDividerMultiplier, assignedRecordDomain, assignedRecordRange);
		int chunkProcessed = perCoreStaticAllocatedChunk;

		clock_t currTS_p8	= clock();
		for(int i=0; i<perCoreStaticAllocatedChunk; i++){
			int selectedChunk = (coreRank-1)*perCoreStaticAllocatedChunk + i;
			doRequestedRecordComp(selectedChunk, bm);
		}

		int selectedChunk;
		while(selectedChunk != -1){
			MPI_Send(&coreRank, 1, MPI_INT, ROOT_ID, WORK_COMPLETE, MPI_COMM_WORLD);
			MPI_Recv(&selectedChunk, 1, MPI_INT, ROOT_ID, WORK_ASSINGED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			if(selectedChunk != -1){
				doRequestedRecordComp(selectedChunk, bm);
				chunkProcessed++;
			}
		}
	}

	if(coreRank == ROOT_ID){
		int received=-1;
		int assignedChunk = (numCores-1)*perCoreStaticAllocatedChunk;
		MPI_Request requests[numCores-1];
		for(int i = 0; i< numCores-1; i++){
			MPI_Irecv(&received, 1, MPI_INT, i+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[i]);
		}
		while(assignedChunk < (numCores-1)*compDividerMultiplier){
			int index_count;
            int indices[numCores-1];

            MPI_Waitsome(numCores-1, requests, &index_count, indices, MPI_STATUSES_IGNORE);

            for(int i = 0; i < index_count; i++)
            {
				if(assignedChunk < (numCores-1)*compDividerMultiplier){
					MPI_Send(&assignedChunk, 1, MPI_INT, indices[i]+1, WORK_ASSINGED, MPI_COMM_WORLD);
					MPI_Irecv(&received, 1, MPI_INT, indices[i]+1, WORK_COMPLETE, MPI_COMM_WORLD, &requests[indices[i]]);
					assignedChunk++;
				}
            }
		}
        MPI_Waitall(numCores-1, requests, MPI_STATUSES_IGNORE);
		int endFlag = -1;
		for(int i=1; i<numCores; i++ ){
			MPI_Send(&endFlag, 1, MPI_INT, i, WORK_ASSINGED, MPI_COMM_WORLD);
		}
	}

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