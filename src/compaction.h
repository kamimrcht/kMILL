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
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, uint k);  /* remove overlaps that are duplicated in vector left, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
vector<edge> removeNotSinglesInRight(const vector<edge>& vect, uint k);  /* remove overlaps that are duplicated in vector right, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
//~ void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR); /* from sequences spotted, get reads' index for reads we do not want the prefix to be present in the next pass */ 
//~ void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR); /* from sequences spotted, get reads' index for reads we do not want the suffix to be present in the next pass */ 
bool compareEdgeByString(const edge& seqL, const edge& seqR);
string compaction(const readStruct& seq1, const readStruct& seq2, uint k);  /*  compaction of two unitigs if they have a k-overlap */
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k);  /*  Recursive func to compact unitigs (if they are sequences or indexes)). If two unitigs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate. */
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k);  /*  checks from the suffixes and prefixes of pairs of unitigs if they can be compacted. Appends overlaps that should be removed in the next pass. */
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);  /* fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs */
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k);  /* fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs */
void cleanDuplicatesInreadStructs(vector <readStruct>& vec);  /* remove duplicates in reads*/
void setreadStructsIndex(vector <readStruct>& vec);  /* give a proper index to reads */
void initVectofreadStructs(vector <readStruct>& vec, const string& sequence);


#endif
