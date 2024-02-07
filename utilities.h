#ifndef UTILITIES_H
#define UTILITIES_H

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <ctime> 

using namespace std;

void getFormattedDataFromCSV(string& file_path, vector<vector<string> > &vec2D, int &attributes, int &totalRecords) {
    string line;
    ifstream records(file_path);
    int ind = 0;

    while (getline (records, line)) {
        vector<string> result;
        boost::split(result, line, boost::is_any_of(","));
        vector<string> vec;
		for(int i=0; i<result.size(); i++) {
			auto last = std::remove_if(result[i].begin(), result[i].end(), [](auto ch) {
        								return ::ispunct(ch) || ::iswpunct(ch);
    								});
			result[i].erase(last, result[i].end()); //Remove junk left by remove_if() at the end of iterator
			boost::to_lower(result[i]);
			vec.push_back(result[i]);
		}
        vec2D.push_back(vec);
    }
    records.close();

	attributes = vec2D[0].size();
	totalRecords = vec2D.size();
}


void getCombinedData(vector<vector<string> > &vec2D, vector<pair<int, string> > &combinedData, int &lenMax) {
	string strSample(50, '0');
	combinedData.resize(vec2D.size());
    int attributes = vec2D[0].size();
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
        for(int j=1; j<attributes; j++) {
            p.second += vec2D[i][j];
        }
		combinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;
	// Padding to make all characters same size
    for (int i = 0; i < vec2D.size(); ++i) {
		int lenDiff		= lenMax - combinedData[i].second.length();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}


void getExactMatches(int &totalUniqueRecords, vector<vector<int> > &exactmatches, vector<pair<int, string> > &combinedData) {
	vector<int> tempVec;
	tempVec.push_back(combinedData[0].first);
	for (int i = 1; i < combinedData.size(); ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(combinedData[i].first);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(combinedData[i].first);
		}
	}
	exactmatches.push_back(tempVec);
	totalUniqueRecords = exactmatches.size();
	cout << "total exact clusters: " << totalUniqueRecords << endl;
}


void getExactAttrs(vector<vector<int> > &exactmatches, vector<pair<int, string> > &combinedData) {
	vector<int> tempVec;
	tempVec.push_back(combinedData[0].first);
	for (int i = 1; i < combinedData.size(); ++i) {
		if(combinedData[i].second.compare(combinedData[i - 1].second) == 0)
			tempVec.push_back(combinedData[i].first);
		else {
			exactmatches.push_back(tempVec);
			tempVec.clear();
			tempVec.push_back(combinedData[i].first);
		}
	}
	exactmatches.push_back(tempVec);
}


void getCombinedData_attr(vector<vector<string> > &vec2D, vector<pair<int, string> > &combinedData, int &lenMax, int attrIndex) {
	string strSample(50, '0');
	combinedData.resize(vec2D.size());
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
		p.second = vec2D[i][attrIndex];
		combinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;
	// Padding to make all characters same size
    for (int i = 0; i < combinedData.size(); ++i) {
		int lenDiff		= lenMax - combinedData[i].second.length();
		if(lenDiff > 0)
			combinedData[i].second	+= strSample.substr(0, lenDiff);
	}
}


void radixSort_std(int &lenMax, vector<pair<int, string> > &strDataArr){
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


void radixSort_utility(int &lenMax, vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	vector<pair<int, string> > tempArr(numRecords);
	
	for (int i = lenMax - 1; i >= 0; --i) {
		vector<int> countArr(42, 0);
		
		for (int j = 0; j < numRecords; ++j) {
			countArr[(strDataArr[j].second)[i] - 48]++;
		}
		
		for (int k = 1; k < 42; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = numRecords - 1; j >= 0; --j)
			tempArr[--countArr[(strDataArr[j].second)[i] - 48]]	= strDataArr[j];
		
		for (int j = 0; j < numRecords; ++j)
			strDataArr[j] = tempArr[j];
	}
}


void getCombinedData_noPad(vector<vector<string> > &vec2D, vector<pair<int, string> > &combinedData, int &lenMax) {
	combinedData.resize(vec2D.size());
    int attributes = vec2D[0].size();
	int max = 0;
	for (int i = 0; i< vec2D.size(); i++) {
		pair<int, string> p;
		p.first = i;
        for(int j=1; j<attributes; j++) {
            p.second += vec2D[i][j];
        }
		combinedData[i]=p;
		if (max<p.second.size()) {
			max = p.second.size();
		}
	}
	lenMax = max;
}


void radixSort_utility_uneven(int lenMax, vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	// first sort them based on sizes
	// then sort all strings sharing same size
	// place them one after another
	// write a function that sorts strings within a certain region in input vector

	vector<pair<int, string> > tempArr(numRecords);
	vector<pair<int, int> > lenRegions;

	// Step 1:
	int start = 0;
	int end = 0;
	vector<int> countLenArr(lenMax, 0);
	for (int i = 0; i < numRecords; ++i) {
		countLenArr[strDataArr[i].second.size()]++;
	}
	for (int k = 1; k < lenMax; ++k) {
		start = countLenArr[k - 1];
		countLenArr[k]	+= countLenArr[k - 1];
		end = countLenArr[k];
		if(start != end) {
			pair<int, int> p;
			p.first = start;
			p.second = end;
			lenRegions.push_back(p);
		}
	}
	for (int i = numRecords - 1; i >= 0; --i) {
			tempArr[--countLenArr[strDataArr[i].second.size()]]	= strDataArr[i];
	}

	for (int i = 0; i < numRecords; ++i){
		strDataArr[i] = tempArr[i];
	}

	// Next Step: Do sorting for each segments

	for(int s=0; s<lenRegions.size(); s++){
		int start = lenRegions[s].first;
		int end = lenRegions[s].second;
		int curLen = strDataArr[start].second.size();
		for (int i = curLen - 1; i >= 0; --i) {
			vector<int> countArr(128, 0);
			countArr[0] = start;
			for (int j = start; j < end; ++j) {
				countArr[(strDataArr[j].second)[i]]++;
			}

			for (int k = 1; k < 128; ++k)
				countArr[k]	+= countArr[k - 1];
			
			for (int j = end - 1; j >= start; --j)
				tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
			
			for (int j = start; j < end; ++j)
				strDataArr[j] = tempArr[j];
		}
	}
}


void radixSort_utility_uneven_improved(int lenMax, vector<pair<int, string> > &strDataArr){
	int numRecords = strDataArr.size();
	// first sort them based on sizes
	// then sort all strings sharing same size
	// place them one after another
	// write a function that sorts strings within a certain region in input vector
	clock_t currTS_p1	= clock();
	vector<pair<int, string> > tempArr(numRecords);
	vector<pair<int, int> > lenRegions;
	vector<int> regionLens;

	// Step 1:
	int start = 0;
	int end = 0;
	vector<int> countLenArr(lenMax+1, 0);
	for (int i = 0; i < numRecords; ++i) {
		countLenArr[strDataArr[i].second.size()]++;
	}
	for (int k = 1; k <= lenMax; ++k) {
		start = countLenArr[k - 1];
		countLenArr[k]	+= countLenArr[k - 1];
		end = countLenArr[k];
		if(start != end) {
			pair<int, int> p;
			p.first = start;
			p.second = end;
			lenRegions.push_back(p);
			regionLens.push_back(k);
		}
	}
	for (int i = numRecords - 1; i >= 0; --i) {
			tempArr[--countLenArr[strDataArr[i].second.size()]]	= strDataArr[i];
	}

	for (int i = 0; i < numRecords; ++i){
		strDataArr[i] = tempArr[i];
	}
	cout<< "Lenmax: " << lenMax << endl;
	double step_1	= (double)(clock() - currTS_p1) / CLOCKS_PER_SEC;
	cout<< "Step 1 Time "<< step_1 << endl;
	clock_t currTS_p2	= clock();
	// Next Step: Do sorting for each segments
	int startRegion = regionLens.size()-1;
	// reuse variable
	start = -1;
	for (int i = lenMax - 1; i >= 0; --i) {
		while((regionLens[startRegion]>i) && (startRegion >= 0)){
			startRegion--;
		}
		if(startRegion >= 0) {
			start = lenRegions[startRegion].second;
		} else {
			start = 0;
		}
		
		// cout<<"At char: "<< i << " start Reg: " << startRegion << "Region lens: " << regionLens[startRegion]<< endl;
		vector<int> countArr(128, 0);
		countArr[0] = start;
		for (int j = start; j < numRecords; ++j) {
			countArr[(strDataArr[j].second)[i]]++;
		}

		for (int k = 1; k < 128; ++k)
			countArr[k]	+= countArr[k - 1];
		
		for (int j = numRecords - 1; j >= start; --j)
			tempArr[--countArr[(strDataArr[j].second)[i]]]	= strDataArr[j];
		
		for (int j = start; j < numRecords; ++j)
			strDataArr[j] = tempArr[j];
	}
	// for(int i=0; i<numRecords-1; i++){
	// 	if(strDataArr[i].second > strDataArr[i+1].second){
	// 		cout<< "Error: " << strDataArr[i].second << " and " << strDataArr[i+1].second << endl;
	// 	}
	// }
	double step_2	= (double)(clock() - currTS_p2) / CLOCKS_PER_SEC;
	cout<< "Step 2 Time "<< step_2 << endl;
}


#endif

