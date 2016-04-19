#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "compaction.h"

#ifndef UTILS
#define UTILS

using namespace std;

bool isCanonical(const string& seq);
char revCompChar(char c);
string revComp(const string& seq);
string getCanonical(const string& seq);
/* read generation */
void createinputlm(uint lr,uint k);
void perfectsReadsFromRef(const string& refName,uint length,uint nbRead);
void randGenome(uint size);
/* debug */
bool isSubSequenceInSequence(const string& subseq, const string& seq);
void sequences2dot(vector<readStruct>& seqV, uint k, unordered_set<uint>& colorNodePref, unordered_set<uint>& colorNodeSuff, unordered_map<uint, uint>& sizesNode);
#endif
