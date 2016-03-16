#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>

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


struct compareEdge{
    bool operator()(const edge& seqL, const edge& seqR){
        return seqL.sequence <seqR.sequence;
    }
};


struct compareRead{
    bool operator()(const readStruct& seqL, const readStruct& seqR){
        return seqL.sequence <seqR.sequence;
    }
};


edge nPrefix(uint n, uint index, const string& sequence);
edge nSuffix(uint n, uint index, const string& sequence);
//~ vector<edge> removeNotSingles(const vector<edge>& vect, unordered_set<string>& edgesToRemove);
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& edgesToRemove, uint k);
vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& edgesToRemove, uint k);
vector<edge> removeEdgesForNextK(const vector<edge>& vect, const unordered_set<string>& edgesToRemove);
bool compareEdgeByString(const edge& seqL, const edge& seqR);
string compaction(const readStruct& seq1, const readStruct& seq2, uint k);
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k);
//~ void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k);
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& edgesToRemove);
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);
void cleanDuplicatesInreadStructs(vector <readStruct>& vec);
void setreadStructsIndex(vector <readStruct>& vec);
void initVectofreadStructs(vector <readStruct>& vec, const string& sequence);


#endif
