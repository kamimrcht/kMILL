#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

#ifndef UTILS
#define UTILS

using namespace std;


string revComp(const string& seq);
string getCanonical(const string& seq);
void createinputlm(uint lr,uint k);
void prefectsReadsFromRef(const string& refName,uint length,uint nbRead);

#endif
