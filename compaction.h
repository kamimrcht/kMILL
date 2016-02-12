#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

#ifndef COMPACTION
#define COMPACTION

using namespace std;

struct edge{
  uint index;
  string sequence;
};


struct readStruct{
  uint index;
  string sequence;
};


struct compareRead{
    bool operator()(const readStruct& seqL, const readStruct& seqR){
        return seqL.sequence <seqR.sequence;
    }
};

edge nPrefix(uint n, uint index, const string& sequence);
edge nSuffix(uint n, uint index, const string& sequence);
vector<edge> removeDuplicates(const vector<edge>& vect);
bool compareEdgeByString(const edge& seqL, const edge& seqR);
string compaction(const readStruct& seq1, const readStruct& seq2, uint k);
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k);
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, int k);
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, int k);
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, int k);
void cleanDuplicatesInreadStructs(vector <readStruct>& vec);
void setreadStructsIndex(vector <readStruct>& vec);
void initVectofreadStructs(vector <readStruct>& vec, string sequence);


#endif
