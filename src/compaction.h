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
  bool takeSuff;
  bool takePref;
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



struct compareRead2{
    bool operator()(const readStruct& seqL, const readStruct& seqR){
        return seqL.sequence <seqR.sequence;
    }
};



edge nPrefix(uint n, uint index, const string& sequence, bool canon);
edge nSuffix(uint n, uint index, const string& sequence, bool canon);
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff);  /* remove overlaps that are duplicated in vector left, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff);  /* remove overlaps that are duplicated in vector right, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove); /* from sequences spotted, get reads' index for reads we do not want the prefix to be present in the next pass */
void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove); /* from sequences spotted, get reads' index for reads we do not want the suffix to be present in the next pass */
bool compareEdgeByString(const edge& seqL, const edge& seqR);
string compaction( readStruct& seq1, const readStruct& seq2, uint k);  /*  compaction of two unitigs if they have a k-overlap */
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k);  /*  Recursive func to compact unitigs (if they are sequences or indexes)). If two unitigs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate. */
void parseVector(vector<vector<edge>> & left, vector<vector<edge>>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemoveInSuff, unordered_set<string>& seqsToRemoveInPref);  /*  checks from the suffixes and prefixes of pairs of unitigs if they can be compacted. Appends overlaps that should be removed in the next pass. */
void fillPrefVector(vector<vector <edge>>& vecLeft, vector<vector <edge>>& vecRight, const readStruct& seq, uint k, string& rev, string& canonPrefix);  /* fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs */
void fillSuffVector(vector<vector <edge>>& vecLeft, vector<vector <edge>>& vecRight, const readStruct& seq, uint k);  /* fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs */
void cleanDuplicatesInreadStructs(vector <readStruct>& vec);  /* remove duplicates in reads*/
void cleanDuplicatesInreadStructs2(vector <readStruct>& vec,uint lol);  /* remove duplicates in reads*/
void setreadStructsIndex(vector <readStruct>& vec);  /* give a proper index to reads */
void initVectofreadStructs(vector <readStruct>& vec, const string& sequence);
void printReadStructsIndex(vector <readStruct>& vec,const string& outfileName);
void fillSuffVector(vector<vector <edge>>& vecLeft, vector<vector <edge>>& vecRight, const readStruct& seq, uint k,string& rev, string& canon);
void fillPrefVector(vector<vector <edge>>& vecLeft, vector<vector <edge>>& vecRight, const readStruct& seq, uint k,string& rev, string& canonPrefix);



#endif
