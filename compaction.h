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
  bool canonical;
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


edge nPrefix(uint n, uint index, const string& sequence, bool canon);
edge nSuffix(uint n, uint index, const string& sequence, bool canon);
//~ vector<edge> removeNotSingles(const vector<edge>& vect, unordered_set<string>& edgesToRemove);
//~ vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& edgesToRemove, uint k);
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff, uint k);
//~ vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& edgesToRemove, uint k);
vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff, uint k);
//~ vector<edge> removeEdgesForNextK(const vector<edge>& vect, const unordered_set<string>& edgesToRemove);
//~ void appendListReadsToRemove(const vector<edge>& vect, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove);
void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove);
void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove);
bool compareEdgeByString(const edge& seqL, const edge& seqR);
//~ string compaction(const readStruct& seq1, const readStruct& seq2, uint k);
string compaction(const readStruct& seq1, const readStruct& seq2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff);
//~ void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k);
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff);

//~ void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k);
//~ void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& edgesToRemove);
//~ void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove);
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemoveInSuff, unordered_set<string>& seqsToRemoveInPref, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff);
//~ void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemovePref);
//~ void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemoveSuff);
void cleanDuplicatesInreadStructs(vector <readStruct>& vec);
void setreadStructsIndex(vector <readStruct>& vec);
void initVectofreadStructs(vector <readStruct>& vec, const string& sequence);


#endif
