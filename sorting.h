#ifndef SORTING_H
#define SORTING_H

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <unordered_set>
#include <queue>
#include <mpi.h>
#include <pthread.h>

using namespace std;

#define SEND_INDICES_TOBE_SORTED 101
#define SEND_STRINGS_TOBE_SORTED 102
#define SEND_SORTED_INDICES 103


void radixSort(int &lenMax, vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string> > tempArr(numRecords);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(256, 0);
		
		for (int j = 0; j < numRecords; ++j) {
			countArr[(strDataArr[j].second)[i]]++;
		}
		
		for (int k = 1; k < 256; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = numRecords - 1; j >= 0; --j)
			tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
		
		for (int j = 0; j < numRecords; ++j)
			strDataArr[j] = tempArr[j];
	}
}

// SampleSort Functions

int getSampleSize(int n, int p){
	double res = log(n)/log(10);
	double s = res;
	// cout<< "Res: " << s << endl;
	while (s < p){
		s = 2*res;
	}
	return (int)s;
}


void generateRandomSample(int s, int* sampleArr, int totalRecords){ 
    std::random_device rd;  // a seed source for the random number engine
    std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> u_distrib(0, totalRecords-1);
 
	int ind = 0;
	std::unordered_set<int> elems(s);
	for (size_t r = 0; r < s;) {
		int v = u_distrib(gen);
		if (auto iter = elems.find(v); iter == elems.end()){
			elems.insert(v);
			sampleArr[ind] = v;
			r++;
			ind++;
		}
	}
}


void getRandomSamples(vector<pair<int, string>> &sampledRecords, int* sampleArr, vector<pair<int, string>> &recordVector, int sampleSize){
	sampledRecords.resize(sampleSize);
	for(int i = 0; i<sampleSize; i++) {
		sampledRecords[i].first = recordVector[sampleArr[i]].first;
		sampledRecords[i].second = recordVector[sampleArr[i]].second;
	}
}


void sortRandomSamples(vector<pair<int, string>> &samples, int* sampleArr, int lenMax){
	radixSort(lenMax, samples);
	// for(int i=0; i<samples.size(); i++){
	// 	sampleArr[i] = samples[i].first;
	// }
}

// does a binary search in sorted samples and returns bucket id
int findBucket(string &target, vector<pair<int, string>> &sortedSamples){
	int l = 0;
	int r = sortedSamples.size()-1;
	while (l<=r){
		int m = floor((l+r)/2);
		int val = target.compare(sortedSamples[m].second);
		if(val > 0){
			l = m+1;
		} else if (val < 0){
			r = m-1;
		} else {
			return m;
		}
	}
    return l;
}


void doBinarySearchAndFindBuckets(vector<pair<int,string>> &sortedSamples, vector<pair<int, string> > &strDataArr, vector<vector<int>> &bucketInds, vector<int> &bucketStarts){
	bucketInds.resize(sortedSamples.size()+1);
	for(int i=0; i<strDataArr.size(); i++){
		int regionInd = findBucket(strDataArr[i].second, sortedSamples);
		// if(regionInd<0 || regionInd > sortedSamples.size()) {
		// 	cout<< "Bucket error " << regionInd << endl;
		// 	break;
		// }
		bucketInds[regionInd].push_back(strDataArr[i].first);
	}
	int sum = 0;

	bucketStarts.resize(bucketInds.size());
	for(int i=0; i<sortedSamples.size()+1; i++){
		bucketStarts[i] = sum;
		sum += bucketInds[i].size();
	}
}

void getBucketAssignment(vector<vector<int>> &bucketInds, vector<vector<int>> &bucketAssignments, vector<int> &totalRecsAssigned, int sampleSize){
	bucketAssignments.resize(sampleSize+1);
	priority_queue<pair<int, int>> pq;
	for(int i=0; i<sampleSize+1; i++){
		pair<int, int> p;
		p.first = INT16_MAX;
		p.second = i;
		pq.push(p);
	}
	for(int i=0; i<bucketInds.size(); i++){
		pair<int, int> p = pq.top();
		pq.pop();
		bucketAssignments[p.second].push_back(i);
		p.first = p.first - bucketInds[i].size();
		pq.push(p);
	}

	totalRecsAssigned.resize(bucketAssignments.size());
	for(int i=0; i< bucketAssignments.size(); i++){
		// cout<<"Buckets for: " << i << endl;
		int sum = 0;
		for(int j=0; j<bucketAssignments[i].size(); j++){
			// cout<< " " << bucketAssignments[i][j] << " ";
			sum += bucketInds[bucketAssignments[i][j]].size();
		}
		totalRecsAssigned[i] = sum;
		// cout<<"Sum: " << sum << endl;
	}
}

void doSorting(vector<vector<int>> &bucketInds, vector<vector<int>> &bucketAssignments, vector<pair<int, string> > &strDataArr, vector<int> &totalRecsAssigned, vector<int> &bucketStarts, int lenMax){
	vector<pair<int, string> > sortedStrings;
	sortedStrings.resize(strDataArr.size());
	for(int i=0; i<bucketAssignments.size(); i++){
		// Sort Step
		vector<pair<int, string> > tempStrings;
		tempStrings.resize(totalRecsAssigned[i]);
		int ind = 0;
		for(int j=0; j<bucketAssignments[i].size(); j++){
			for(int k = 0; k< bucketInds[bucketAssignments[i][j]].size(); k++){
				tempStrings[ind] = strDataArr[bucketInds[bucketAssignments[i][j]][k]];
				ind++;
			}
		}
		radixSort(lenMax, tempStrings);
		// for(int u =1; u<tempStrings.size(); u++){
		// 	if (tempStrings[u].second.compare(tempStrings[u-1].second)<0){
		// 		cout<< "ERROR" << endl;
		// 		cout<< tempStrings[u-1].second << endl;;
		// 		cout<< tempStrings[u].second << endl;
		// 	}
		// }

		// Merge Step
		ind = 0;
		for(int j=0; j<bucketAssignments[i].size(); j++){
			int ind_start = bucketStarts[bucketAssignments[i][j]];
			for(int k =0; k<bucketInds[bucketAssignments[i][j]].size(); k++){
				sortedStrings[ind_start+k] = tempStrings[ind];
				// cout<< tempStrings[ind].second << endl;
				ind++;
			}
		}
	}
	// for(int u =1; u<sortedStrings.size(); u++){
	// 	if (sortedStrings[u].second.compare(sortedStrings[u-1].second)<0){
	// 		cout<< "ERROR" << endl;
	// 		cout<< sortedStrings[u-1].second << endl;;
	// 		cout<< sortedStrings[u].second << endl;
	// 	}
	// }
	strDataArr = sortedStrings;

}

// void parallelSampleSort_1(int lenMax, int numCores, vector<pair<int, string> > &strDataArr, int coreRank){
// 	vector<pair<int, string>> sortedSampledRecords;
// 	vector<pair<int, string>> chosenSortedSampledRecords;
// 	vector<vector<int>> bucketInds;
// 	vector<vector<int>> bucketAssignments;
// 	vector<int> bucketStarts;
// 	vector<int> totalRecsAssigned;
// 	int sampleSize;
// 	int lenSize;
// 	// int factor = 40;
// 	// int sampleSize = getSampleSize(strDataArr.size(), numCores-1);

// 	if(coreRank == 0){
// 		sampleSize = 2*(numCores-1);
// 		lenSize = lenMax;
// 	}

// 	int* tRecs = new int[numCores];

// 	MPI_Bcast(&lenSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
// 	MPI_Bcast(&sampleSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


// 	if (coreRank == 0) {
// 		cout<< "SampleSize: " << sampleSize << endl;
// 		int regionPerSample = 10;
// 		int effectiveSampleSize = (sampleSize + 1)*regionPerSample;
// 		int* sampleInd;
// 		sampleInd = new int[effectiveSampleSize];
// 		generateRandomSample(effectiveSampleSize, sampleInd, strDataArr.size());
// 		getRandomSamples(sortedSampledRecords, sampleInd, strDataArr, effectiveSampleSize);
// 		sortRandomSamples(sortedSampledRecords, sampleInd, lenSize);
// 		for(int i=regionPerSample-1; i< effectiveSampleSize-regionPerSample+1; i++){
// 			if((i%(regionPerSample-1)) == 0){
// 				chosenSortedSampledRecords.push_back(sortedSampledRecords[i]);
// 			}
// 		}
// 		doBinarySearchAndFindBuckets(chosenSortedSampledRecords, strDataArr, bucketInds, bucketStarts);
// 		getBucketAssignment(bucketInds, bucketAssignments, totalRecsAssigned, numCores);
// 		for(int i=0; i<numCores-1; i++){
// 			tRecs[i] = totalRecsAssigned[i];
// 		}
// 	}

// 	MPI_Bcast(tRecs, numCores-1, MPI_INT, 0, MPI_COMM_WORLD);


// 	if(coreRank == 0){
// 		for(int i=0; i<numCores-1; i++){
// 			int msgCharSize = totalRecsAssigned[i]*(lenMax+1);
// 			char *sortSuperString = new char[msgCharSize];
// 			int *recIds = new int[totalRecsAssigned[i]];
// 			int ind = 0;
// 			int recIND = 0;
// 			int totalTokens = 0;
// 			for(int j=0; j<bucketAssignments[i].size(); j++){
// 				for(int k = 0; k< bucketInds[bucketAssignments[i][j]].size(); k++){
// 					for(size_t l=0; l< std::strlen(strDataArr[bucketInds[bucketAssignments[i][j]][k]].second.c_str()); l++){
// 						sortSuperString[ind++] = strDataArr[bucketInds[bucketAssignments[i][j]][k]].second[l];
// 					}
// 					sortSuperString[ind++] = '|';
// 					totalTokens++;
// 					recIds[recIND++] = bucketInds[bucketAssignments[i][j]][k];
// 				}
// 			}
// 			// cout<< "Total Records Assigned to core "<< i+1 << " is " << totalRecsAssigned[i] << endl;
// 			MPI_Send(recIds, totalRecsAssigned[i], MPI_INT, i+1, SEND_INDICES_TOBE_SORTED, MPI_COMM_WORLD);
// 			MPI_Send(sortSuperString, msgCharSize, MPI_CHAR, i+1, SEND_STRINGS_TOBE_SORTED, MPI_COMM_WORLD);
// 		}
// 		// cout<< "DATA SENT FOR SORTING" << endl;
// 	}

// 	if(coreRank != 0) {
// 		cout<< "CoreRank" << coreRank << " Has " << tRecs[coreRank-1] << " Records" << endl;
// 		int msgCharSize = tRecs[coreRank-1]*(lenSize+1);
// 		char *sortSuperString = new char[msgCharSize];
// 		int *recIds = new int[tRecs[coreRank-1]];

// 		MPI_Recv(recIds, tRecs[coreRank-1], MPI_INT, 0, SEND_INDICES_TOBE_SORTED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
// 		MPI_Recv(sortSuperString, msgCharSize, MPI_CHAR, 0, SEND_STRINGS_TOBE_SORTED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

// 		// cout<< "CoreRank " << coreRank << " RECEIVED ALL DATA" << endl;
// 		// for(int i=0; i<10; i++){
// 		// 	cout<< "Printing recIDS " << recIds[i] << endl; 
// 		// }
// 		string superStringStr(sortSuperString);

// 		vector<pair<int, string> > sortedStrings;
// 		sortedStrings.resize(tRecs[coreRank-1]);
	
// 		vector<string> attrStrs;
// 		boost::split(attrStrs, superStringStr, boost::is_any_of("|"));
// 		for(int i=0; i< tRecs[coreRank-1]; i++){
// 			pair<int, string> p;
// 			p.first = recIds[i];
// 			p.second = attrStrs[i];
// 			sortedStrings[i] = p;
// 		}

// 		radixSort(lenSize, sortedStrings);

// 		for(int i=0; i<tRecs[coreRank-1]; i++){
// 			recIds[i] = sortedStrings[i].first;
// 		}
// 		// cout<< "READY TO SEND OUT SORTED STRINGS" << endl;
// 		MPI_Send(recIds, tRecs[coreRank-1], MPI_INT, 0, SEND_SORTED_INDICES, MPI_COMM_WORLD);
// 	}
	
// 	if(coreRank == 0){
// 		vector<pair<int, string> > sortedStrings;
// 		sortedStrings.resize(strDataArr.size());

// 		for(int i=0; i<numCores-1; i++){
// 			int *recIds = new int[totalRecsAssigned[i]];
// 			MPI_Recv(recIds, totalRecsAssigned[i], MPI_INT, i+1, SEND_SORTED_INDICES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

// 			int ind = 0;
// 			for(int j=0; j<bucketAssignments[i].size(); j++){
// 				int ind_start = bucketStarts[bucketAssignments[i][j]];
// 				for(int k =0; k<bucketInds[bucketAssignments[i][j]].size(); k++){
// 					// cout<< "Sorted Indice: " << recIds[ind] << endl;
// 					// cout<< "Sorted String: " << strDataArr[recIds[ind]].second << endl;
// 					sortedStrings[ind_start+k] = strDataArr[recIds[ind]];
// 					ind++;
// 				}
// 			}
// 		}
// 		strDataArr = sortedStrings;
// 	}
// }

void parallelSampleSort_1(int lenMax, int numCores, vector<pair<int, string> > &strDataArr, int coreRank){
	vector<pair<int, string>> sortedSampledRecords;
	vector<pair<int, string>> chosenSortedSampledRecords;
	vector<vector<int>> bucketInds;
	vector<vector<int>> bucketAssignments;
	vector<int> bucketStarts;
	vector<int> totalRecsAssigned;
	int sampleSize;
	int lenSize;
	// int factor = 40;
	// int sampleSize = getSampleSize(strDataArr.size(), numCores-1);

	if(coreRank == 0){
		sampleSize = 2*(numCores-1);
		lenSize = lenMax;
	}

	int* tRecs = new int[numCores];

	MPI_Bcast(&lenSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&sampleSize, 1, MPI_INT, 0, MPI_COMM_WORLD);


	if (coreRank == 0) {
		cout<< "SampleSize: " << sampleSize << endl;
		int regionPerSample = 10;
		int effectiveSampleSize = (sampleSize + 1)*regionPerSample;
		int* sampleInd;
		sampleInd = new int[effectiveSampleSize];
		generateRandomSample(effectiveSampleSize, sampleInd, strDataArr.size());
		getRandomSamples(sortedSampledRecords, sampleInd, strDataArr, effectiveSampleSize);
		sortRandomSamples(sortedSampledRecords, sampleInd, lenSize);
		for(int i=regionPerSample-1; i< effectiveSampleSize-regionPerSample+1; i++){
			if((i%(regionPerSample-1)) == 0){
				chosenSortedSampledRecords.push_back(sortedSampledRecords[i]);
			}
		}
		doBinarySearchAndFindBuckets(chosenSortedSampledRecords, strDataArr, bucketInds, bucketStarts);
		getBucketAssignment(bucketInds, bucketAssignments, totalRecsAssigned, numCores);
		for(int i=0; i<numCores-1; i++){
			tRecs[i] = totalRecsAssigned[i];
		}
	}

	MPI_Bcast(tRecs, numCores-1, MPI_INT, 0, MPI_COMM_WORLD);


	if(coreRank == 0){
		for(int i=0; i<numCores-1; i++){
			int msgCharSize = totalRecsAssigned[i]*(lenMax+1);
			char *sortSuperString = new char[msgCharSize];
			int *recIds = new int[totalRecsAssigned[i]];
			int ind = 0;
			int recIND = 0;
			int totalTokens = 0;
			for(int j=0; j<bucketAssignments[i].size(); j++){
				for(int k = 0; k< bucketInds[bucketAssignments[i][j]].size(); k++){
					for(size_t l=0; l< std::strlen(strDataArr[bucketInds[bucketAssignments[i][j]][k]].second.c_str()); l++){
						sortSuperString[ind++] = strDataArr[bucketInds[bucketAssignments[i][j]][k]].second[l];
					}
					sortSuperString[ind++] = '|';
					totalTokens++;
					recIds[recIND++] = bucketInds[bucketAssignments[i][j]][k];
				}
			}
			// cout<< "Total Records Assigned to core "<< i+1 << " is " << totalRecsAssigned[i] << endl;
			MPI_Send(recIds, totalRecsAssigned[i], MPI_INT, i+1, SEND_INDICES_TOBE_SORTED, MPI_COMM_WORLD);
			MPI_Send(sortSuperString, msgCharSize, MPI_CHAR, i+1, SEND_STRINGS_TOBE_SORTED, MPI_COMM_WORLD);
		}
		// cout<< "DATA SENT FOR SORTING" << endl;
	}

	if(coreRank != 0) {
		cout<< "CoreRank" << coreRank << " Has " << tRecs[coreRank-1] << " Records" << endl;
		int msgCharSize = tRecs[coreRank-1]*(lenSize+1);
		char *sortSuperString = new char[msgCharSize];
		int *recIds = new int[tRecs[coreRank-1]];

		MPI_Recv(recIds, tRecs[coreRank-1], MPI_INT, 0, SEND_INDICES_TOBE_SORTED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(sortSuperString, msgCharSize, MPI_CHAR, 0, SEND_STRINGS_TOBE_SORTED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		// cout<< "CoreRank " << coreRank << " RECEIVED ALL DATA" << endl;
		// for(int i=0; i<10; i++){
		// 	cout<< "Printing recIDS " << recIds[i] << endl; 
		// }
		string superStringStr(sortSuperString);

		vector<pair<int, string> > sortedStrings;
		sortedStrings.resize(tRecs[coreRank-1]);
	
		vector<string> attrStrs;
		boost::split(attrStrs, superStringStr, boost::is_any_of("|"));
		for(int i=0; i< tRecs[coreRank-1]; i++){
			pair<int, string> p;
			p.first = recIds[i];
			p.second = attrStrs[i];
			sortedStrings[i] = p;
		}

		radixSort(lenSize, sortedStrings);

		for(int i=0; i<tRecs[coreRank-1]; i++){
			recIds[i] = sortedStrings[i].first;
		}
		// cout<< "READY TO SEND OUT SORTED STRINGS" << endl;
		MPI_Send(recIds, tRecs[coreRank-1], MPI_INT, 0, SEND_SORTED_INDICES, MPI_COMM_WORLD);
	}
	
	if(coreRank == 0){
		vector<pair<int, string> > sortedStrings;
		sortedStrings.resize(strDataArr.size());

		for(int i=0; i<numCores-1; i++){
			int *recIds = new int[totalRecsAssigned[i]];
			MPI_Recv(recIds, totalRecsAssigned[i], MPI_INT, i+1, SEND_SORTED_INDICES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int ind = 0;
			for(int j=0; j<bucketAssignments[i].size(); j++){
				int ind_start = bucketStarts[bucketAssignments[i][j]];
				for(int k =0; k<bucketInds[bucketAssignments[i][j]].size(); k++){
					// cout<< "Sorted Indice: " << recIds[ind] << endl;
					// cout<< "Sorted String: " << strDataArr[recIds[ind]].second << endl;
					sortedStrings[ind_start+k] = strDataArr[recIds[ind]];
					ind++;
				}
			}
		}
		strDataArr = sortedStrings;
	}
}


struct arg_struct {
    int rank;
    int threadID;

};


#endif

