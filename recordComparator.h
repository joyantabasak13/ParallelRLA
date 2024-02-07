#ifndef RECORD_COMPARATOR_H
#define RECORD_COMPARATOR_H

#include "unionFind.h"
#include <vector>
#include <iostream>

using namespace std;

class RecordComparator{
    public:
        int threshold = 99;
        int cumulativeDistanceThreshold = 0;
        vector<int> attrDistThreshold;
        int matArr[50][50] = {0};
        int attributes = 0;

        void setComparator(int gThreshold, vector<int>& perAttrThreshold){
            this->cumulativeDistanceThreshold = gThreshold;
            this->attrDistThreshold = perAttrThreshold;
            this->attributes = attrDistThreshold.size();
        }

        int calculateBasicED2(string &str1, string &str2, int threshRem){
            int row, col, i, j;
            row = str1.length() + 1;
            col = str2.length() + 1;

            for (i = 0; i < row; i++)
            {
                for (j = 0; j < col; j++)
                {
                    if (i == 0)
                        matArr[i][j] = j;
                    else if (j == 0)
                        matArr[i][j] = i;
                    else
                    {
                        if ((int)str1[i - 1] == (int)str2[j - 1])
                            matArr[i][j] = min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
                        else
                            matArr[i][j] = min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
                            // if(i>1 && j>1) {
                            // 	if((str1[i-1]==str2[j-2]) && (str1[i-2] == str2[j-1])) {
                            // 		matArr[i][j] = min(matArr[i][j], matArr[i-2][j-2]+1);
                            // 	}
                            // }
                    }

                    if ((row - col) == (i - j) && (matArr[i][j] > threshRem))
                    {
                        return threshold + 1;
                    }
                }
            }
            
            return (matArr[row - 1][col - 1]);
        }

        int calculateBasicED(string &str1, string &str2, int threshRem){
            int dist = threshRem;

            if (abs((int)(str1.length() - str2.length())) > dist)
                return threshold + 1;
            else if ((2 * dist + 1) >= max(str1.length(), str2.length()))
                return calculateBasicED2(str1, str2, dist);
            else
            {
                
                string s1, s2;
                int row, col, diagonal;
                int i, j;

                if (str1.length() > str2.length()){
                    s1 = str2;
                    s2 = str1;
                } else {
                    s1 = str1;
                    s2 = str2;
                }

                row = s1.length() + 1;
                col = 2 * dist + 1;
                diagonal = dist + s2.length() - s1.length();

                for (i = 0; i < dist + 1; i++)
                {
                    for (j = dist - i; j < col; j++)
                    {
                        if (i == 0)
                            matArr[i][j] = j - dist;
                        else if (j == (dist - i))
                            matArr[i][j] = matArr[i - 1][j + 1] + 1;
                        else if (j != (col - 1))
                        {
                            if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
                                matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                            else {
                                matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                            }
                                    
                        }
                        else
                        {
                            if ((int)s1[i - 1] == (int)s2[j - (dist - i) - 1])
                                matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
                            else
                                matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
                        }

                        if ((j == diagonal) && matArr[i][j] > dist)
                            return threshold + 1;
                    }
                }

                for (i = dist + 1; i < s2.length() - dist + 1; i++)
                {
                    for (j = 0; j < col; j++)
                    {
                        if (j == 0)
                        {
                            if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
                                matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
                            else
                                matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
                        }
                        else if (j != (col - 1))
                        {
                            if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
                                matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                            else
                                matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                        }
                        else
                        {
                            if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
                                matArr[i][j] = min(matArr[i - 1][j], matArr[i][j - 1] + 1);
                            else
                                matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
                        }
                        if ((j == diagonal) && (matArr[i][j] > dist))
                            return threshold + 1;
                    }
                }

                for (i = s2.length() - dist + 1; i < row; i++)
                {
                    for (j = 0; j < col - i + s2.length() - dist; j++)
                    {
                        if (j == 0)
                        {
                            if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
                                matArr[i][j] = min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
                            else
                                matArr[i][j] = min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
                        }
                        else
                        {
                            if ((int)s1[i - 1] == (int)s2[j + (i - dist) - 1])
                                matArr[i][j] = min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                            else
                                matArr[i][j] = min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
                        }
                        if ((j == diagonal) && (matArr[i][j] > dist))
                            return threshold + 1;
                    }
                }
                return matArr[row - 1][diagonal];
            }
        }

        bool isLinkageOk(vector<string> &a, vector<string> &b)
        {
            int cumulativeDist = 0;
            for (int i = 1; i < attributes+1; i++)
            {	
                int dist = calculateBasicED(a[i], b[i], this->attrDistThreshold[i-1]);
                if (dist > this->attrDistThreshold[i-1]){
                    return false;
                } else {
                    cumulativeDist += dist;
                    if (cumulativeDist > this->cumulativeDistanceThreshold){
                        return false;
                    }
                }
            }
            return true;
        }

        bool isLinkageOk(string &a, string &b)
        {
            int dist = calculateBasicED(a, b, this->cumulativeDistanceThreshold);
            if (dist > this->cumulativeDistanceThreshold){
                return false;
            }
            return true;
        }

        bool isLinkageOk(vector<string> &a, vector<string> &b, int rec_i, int rec_j, vector<vector<int>> &parentArrVec)
        {
            int cumulativeDist = 0;
            for (int i = 1; i < parentArrVec.size(); i++)
            {	
                if (parentArrVec[i][rec_i] != parentArrVec[i][rec_j]){
                    return false;
                }
            }
            return this->isLinkageOk(a,b);
        }
};



#endif