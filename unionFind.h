#ifndef UNIONFIND_H
#define UNIONFIND_H

#include <vector>

using namespace std;

class UnionFind {
    public:
        int numSets = 0;
        int* parentInd;

        UnionFind();
        void setVariable(int numRecords);
        int find(int recID);
        void weightedUnion(int recId1, int recId2);
        bool isConnected(int recId1, int recId2);
        int getSetCount();
        void setParent(int recId, int root);
        void setParentArr(int parentArr[]);
        int getRootVal(int recId);
        void setRootVal(int recId, int val);

};

UnionFind::UnionFind() { // Constructor with parameters
    numSets = 0;
}

void UnionFind::setVariable(int numRecords) {
    this->numSets = numRecords;
    this->parentInd = new int[numRecords];
    for (int i = 0; i < numRecords; i++)
    {
        this->parentInd[i] = -1;
    }
}

int UnionFind::find(int recID) {
    int root = recID;
    while (parentInd[root] >= 0){
        root = parentInd[root];
    }
    while (recID != root) { 
        int newParent = parentInd[recID]; 
        parentInd[recID] = root; 
        recID = newParent; 
    }
    return root;
}

void UnionFind::weightedUnion(int recId1, int recId2) {
    int i = find(recId1); 
    int j = find(recId2); 
    if (i == j) {
        return;
    }
    // make smaller root point to larger one
    if (parentInd[i] > parentInd[j]) { 
        parentInd[j] += parentInd[i]; 
        parentInd[i] = j; 
    }
    else { 
        parentInd[i] += parentInd[j];
        parentInd[j] = i;
    }
    numSets--;
}

bool UnionFind::isConnected(int recId1, int recId2) {
    if (find(recId1) == find(recId2)) {
        return true;
    } else {
        return false;
    }
}

int UnionFind::getSetCount() {
    return this->numSets;
}

void UnionFind::setParent(int recId, int root) {
    this->parentInd[recId] = root;
}

int UnionFind::getRootVal(int recId) {
    int root = find(recId);
    return parentInd[root];
}

void UnionFind::setRootVal(int recId, int val) {
    int root = find(recId);
    parentInd[root] = val;
}

void UnionFind::setParentArr(int parentArr[]){
    this->parentInd = parentArr;
}


#endif

