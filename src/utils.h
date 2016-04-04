#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

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
/* debug */
bool isSubSequenceInSequence(const string& subseq, const string& seq);

#endif
