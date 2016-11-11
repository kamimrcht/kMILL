#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include "readAndSortInputFile.h"
#include "utils.h"
#include <unordered_set>
#include <unordered_map>
#include <set>



using namespace std;



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./pathsCleaner reads.fasta threshold" << endl;
		return 0;
	}
	string fileName = argv[1];
	uint thresholdCleaning = stoi(argv[2]);
	ifstream readStructFile(fileName);
	if(not readStructFile){
		cout<<"No such file ..."<<endl;
		return 1;
	}
	uint nbBuckets(1);
	vector <ofstream> outFiles(nbBuckets);
	string sequence,sequence2;
	vector <readStruct> sequencesVec;
	openBuckets(outFiles);
	createReadBuckets(nbBuckets, readStructFile, outFiles);
	fillSortCleanBuckets(nbBuckets, sequencesVec,thresholdCleaning);
	removeReadFiles(nbBuckets);
	setreadStructsIndex(sequencesVec);
	printReadStructsIndex(sequencesVec,"noduplicate.fa");
	return 0;
}
