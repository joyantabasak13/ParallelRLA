#ifndef BLOCKING_H
#define BLOCKING_H

#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <cstdlib>

using namespace std;

class BlockingMethods{
    public:
        int KMER = 3;
        int alphabetSize = 36;
        int blockingAttrIndex = 1;
        long long int blockIDRange = -1;
        long long int totalUniqueRecords;
        long long int totalCompRequired;
        vector<pair<long long int, int> > blockingIDList;
        vector<pair<long long int, int> > boundaryArr;
        vector<vector<string>> uniqueRecords;
        vector<string> uniqueAttrs;

        void setBlocking(int k, int alpha, int blockAttr, vector<vector<string>> &uRecords){
            KMER = k;
            alphabetSize = alpha;
            blockingAttrIndex = blockAttr;
            uniqueRecords = uRecords;
            totalUniqueRecords = uRecords.size();
        }

        void setBlocking(int k, int alpha, vector<string> &uAttrs){
            KMER = k;
            alphabetSize = alpha;
            blockingAttrIndex = -1;
            uniqueAttrs = uAttrs;
            totalUniqueRecords = uAttrs.size();
        }


        void getNormalBlockingIDArray() {
            string strSample(KMER, 'a');
            blockIDRange = pow(alphabetSize, KMER);
            long long int blockID;
            string blockingStr;
            for (int i = 0; i < totalUniqueRecords; i++) {
                blockingStr = uniqueRecords[i][blockingAttrIndex];
                while(blockingStr.size() < KMER) {                    
                    blockingStr += strSample.substr(0, KMER - blockingStr.length());
                }
                int blkstrSize = blockingStr.size();
                for (int j = 0; j < blkstrSize; ++j) {
                    blockID = 0;
                    for (int k = 0; k < KMER; ++k)
                    {
                        int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
                        if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
                            blockID += (blockStrCharCode - 97) * pow(alphabetSize,k);
                        } else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
                            blockID += (blockStrCharCode - 48 + 26) * pow(alphabetSize,k);
                        }
                    }
                    pair<long long int,int> blkRecPair;
                    blkRecPair.first = blockID;
                    blkRecPair.second = i;
                    blockingIDList.push_back(blkRecPair);

                    // if(blockID<0 || blockID>=blockIDRange) {
                    //     cout<< "Invalid BlockID: " << blockID << endl;
                    //     cout<< "BlockingString: " << blockingStr << endl;
                    // }
                }	
            }
        }

        void getNormalBlockingIDArrayAttrs() {
            string strSample(KMER, 'a');
            blockIDRange = pow(alphabetSize, KMER);
            long long int blockID;
            string blockingStr;
            for (int i = 0; i < totalUniqueRecords; i++) {
                blockingStr = uniqueAttrs[i];
                while(blockingStr.size() < KMER) {                    
                    blockingStr += strSample.substr(0, KMER - blockingStr.length());
                }
                int blkstrSize = blockingStr.size();
                for (int j = 0; j < blkstrSize; ++j) {
                    blockID = 0;
                    for (int k = 0; k < KMER; ++k)
                    {
                        int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
                        if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
                            blockID += (blockStrCharCode - 97) * pow(alphabetSize,k);
                        } else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
                            blockID += (blockStrCharCode - 48 + 26) * pow(alphabetSize,k);
                        }
                    }
                    pair<long long int,int> blkRecPair;
                    blkRecPair.first = blockID;
                    blkRecPair.second = i;
                    blockingIDList.push_back(blkRecPair);

                    // if(blockID<0 || blockID>=blockIDRange) {
                    //     cout<< "Invalid BlockID: " << blockID << endl;
                    //     cout<< "BlockingString: " << blockingStr << endl;
                    // }
                }	
            }
        }

        void getSuperBlockingIDArray() {
            blockIDRange = pow(alphabetSize, KMER+1);
            string strSample(KMER, 'a');
            int blockID;
            string blockingStr;
            int perAplhaBlocks = pow(alphabetSize, KMER);
            for (int i = 0; i < totalUniqueRecords; i++) {
                blockingStr = uniqueRecords[i][blockingAttrIndex];
                string temp_str = uniqueRecords[i][blockingAttrIndex];
                while(blockingStr.size() < KMER) {                    
                    blockingStr += strSample.substr(0, KMER - blockingStr.length());
                }
                int blkstrSize = blockingStr.size();
                for (int j = 0; j < blkstrSize; ++j) {
                    blockID = 0;
                    for (int k = 0; k < KMER; ++k)
                    {
                        int blockStrCharCode = (int)blockingStr[(j+k)%blkstrSize];
                        if ((blockStrCharCode >= 97) && (blockStrCharCode <= 122)) {
                            blockID += (blockStrCharCode - 97) * pow(alphabetSize,k);
                        } else if ((blockStrCharCode >= 48) && (blockStrCharCode <= 57)){
                            blockID += (blockStrCharCode - 48 + 26) * pow(alphabetSize,k);
                        }
                    }
                    pair<int,int> blkRecPair;
                    blockID = (blockingStr[0]-97)*perAplhaBlocks + blockID;
                    blkRecPair.first = blockID;
                    blkRecPair.second = i;
                    blockingIDList.push_back(blkRecPair);

                    // if(blockID<0 || blockID>=blockIDRange) {
                    //     cout<< "Invalid BlockID: " << blockID << endl;
                    //     cout<< "BlockingString: " << blockingStr << endl;
                    // }
                }	
            }
        }


        void sortBlockingIDArray() {
            int numRecords = blockingIDList.size();
            vector<pair<long long int, int>> tempArr(numRecords);
            vector<int> countArr(blockIDRange, 0);
            for (int j = 0; j < numRecords; ++j) {
                countArr[blockingIDList[j].first]++;
            }
            // Do prefix sum
            for (long long int k = 1; k < blockIDRange; ++k)
                countArr[k]	+= countArr[k - 1];
            for (int j = numRecords - 1; j >= 0; --j)
                tempArr[--countArr[blockingIDList[j].first]] = blockingIDList[j];
            for (int j = 0; j < numRecords; ++j)
                blockingIDList[j] = tempArr[j];
        }

        void removeRedundentBlockingID() {
            int numRecords = blockingIDList.size();
            vector<pair<long long int, int>> tempArr;
            long long int totalUniqueBlocks = 1;
            tempArr.push_back(blockingIDList[0]);
            for (int i = 1; i<numRecords; i++) {
                if (blockingIDList[i].first != blockingIDList[i-1].first) {
                    totalUniqueBlocks++;
                }
                if ( ! ((blockingIDList[i].first == blockingIDList[i-1].first) && (blockingIDList[i].second == blockingIDList[i-1].second))) {
                    tempArr.push_back(blockingIDList[i]);
                }
            }
            blockingIDList = tempArr;

            // cout<< "Total Unique Blocks: "<< totalUniqueBlocks << endl;
            // cout<< "Total unique block-rec pairs: " << numRecords << endl;
            // cout<< "Removed redundant block-rec pairs: " << numRecords - blockingIDList.size() << endl; 
        }

        void findBlockBoundaries() {
            //just keep starting ind and range.
            // boundaryArr.resize(totalUniqueBlocks);
            int numRecords = blockingIDList.size();
            totalCompRequired = 0;
            long long int startInd = 0;
            int range = 0;
            int curBlockInd = 0;
            for (long long int i = 1; i<numRecords; i++) {
                if (blockingIDList[i].first != blockingIDList[i-1].first) {
                    range = i-startInd;
                    pair<long long int,int> p;
                    p.first = startInd;
                    p.second = range;
                    boundaryArr.push_back(p);
                    totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
                    curBlockInd++;
                    startInd = i;
                }
            }
            // Enter last Block info
            range = numRecords-startInd;
            totalCompRequired = totalCompRequired + ceil((range*(range-1))/2);
            pair<long long int,int> p;
            p.first = startInd;
            p.second = range;
            boundaryArr.push_back(p);
        }

        void findBlockAssignments(int numCores, int numThreads, int compDividerMultiplier, vector<vector<pair<long long int, int> > > &assignedRecordDomain, vector<vector<pair<long long int, int> > > &assignedRecordRange) {
            int totalChunks = numCores*numThreads*compDividerMultiplier;
            assignedRecordDomain.resize(totalChunks);
            assignedRecordRange.resize(totalChunks);

            long long int compThreshold = ceil(totalCompRequired/totalChunks);
            long long int curAssignmentSize = 0;
            int recDomainEndInd = 0;
            int remRecCount = 0;
            int curChunk = 0;

            for (int j = 0; j < boundaryArr.size(); j++)
            {
                int range = boundaryArr[j].second;

                long long int curBlockSize = ceil((range*(range-1))/2);

                // assign all remaining comparisions to the last core
                if(curChunk == totalChunks - 1) {
                    compThreshold = LLONG_MAX;
                }

                if ((curAssignmentSize+curBlockSize) < compThreshold) {
                    assignedRecordDomain[curChunk].push_back(boundaryArr[j]);
                    assignedRecordRange[curChunk].push_back(boundaryArr[j]);
                    curAssignmentSize = curAssignmentSize + curBlockSize;
                } else {
                    remRecCount = boundaryArr[j].second;
                    while(curAssignmentSize + remRecCount < compThreshold){
                        recDomainEndInd++;
                        curAssignmentSize+=remRecCount-1;
                        remRecCount--;
                    }
                    if(recDomainEndInd != 0){
                        pair<long long int,int> dom;
                        dom.first = boundaryArr[j].first;
                        dom.second = recDomainEndInd;
                        assignedRecordDomain[curChunk].push_back(dom);
                        assignedRecordRange[curChunk].push_back(boundaryArr[j]);
                        pair<long long int, int> newBound;
                        newBound.first = boundaryArr[j].first+recDomainEndInd;
                        newBound.second = boundaryArr[j].second-recDomainEndInd;
                        boundaryArr.push_back(newBound);
                        recDomainEndInd = 0;
                        remRecCount = 0;
                    } else {
                        j--;
                    }
                    curChunk++;
                    // cout<< "Assigned Chunk: " << curChunk << endl;
                    curAssignmentSize = 0;
                }
            }
        }

        
};

#endif


