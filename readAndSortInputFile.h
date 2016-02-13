#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "compaction.h"

#ifndef READANDSORT
#define READANDSORT

using namespace std;


uint64_t transformStringToHash(string read);
void createReadBuckets(uint nbBuckets, ifstream& readStructFile, vector <ofstream>& outFiles);
void openBuckets(vector<ofstream>& outFiles);
void fillSortCleanBuckets(uint nbBuckets, vector <readStruct>& sequencesVec);
void removeReadFiles(uint nbBuckets);

#endif
