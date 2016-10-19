#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "compaction.h"

#ifndef READANDSORT
#define READANDSORT

using namespace std;


uint64_t transformStringToHash(const string& read);
void createReadBuckets(uint nbBuckets, ifstream& readStructFile, vector <ofstream>& outFiles);  /* create buckets to store canonical sequences of the reads */
void openBuckets(vector<ofstream>& outFiles);  /* open buckets */
void fillSortCleanBuckets(uint nbBuckets, vector <readStruct>& sequencesVec, bool);  /* get sequences from all buckets, sort them and remove duplicated reads */
void removeReadFiles(uint nbBuckets);  /* remove bucket files */
//test
void createReadBuckets2ndPass(uint nbBuckets, const vector <string>& sequences2ndPass, vector <ofstream>& outFiles);
#endif
