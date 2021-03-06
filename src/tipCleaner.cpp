#include <fstream>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include "compaction.h"
#include "readAndSortInputFile.h"
#include "utils.h"
#include <unordered_set>
#include <unordered_map>
#include <set>



using namespace std;
uint nBucketsOverlapmain(100);



int main(int argc, char ** argv){
	//ARGUMENT PARSING
	if (argc < 4){
		cout << "command line: ./kMILL reads.fasta k minSize" << endl;
		return 0;
	}
	auto startChrono=chrono::system_clock::now();
	string fileName = argv[1];
	uint k(stoi(argv[2]));
	uint minSize(stoi(argv[3]));
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

	//SEQUENCES DUPLICATE CLEANING
	uint maxSize=createReadBuckets(nbBuckets, readStructFile, outFiles);
	string titre("tiped.fa");
	ofstream out(titre);
	if(not out.good()){
		cout<<"probleme outfile"<<endl;
	}
	fillSortCleanBuckets(nbBuckets, sequencesVec,0);
	removeReadFiles(nbBuckets);
	setreadStructsIndex(sequencesVec);
	cout<<sequencesVec.size()<<endl;
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Init took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
	//SEQUENCES TIPING
	unordered_set<string> seqsToRemoveInSuff,seqsToRemoveInPref;
	vector <vector<edge>> right(nBucketsOverlapmain);  // vector of canonical suffixes
	vector <vector<edge>> left(nBucketsOverlapmain); //  vector of canonical prefixes
	//FILLING
	string rev,canon;
	for (uint i(0); i<sequencesVec.size(); ++i){
		if (sequencesVec[i].sequence.size() >= k){
			fillPrefVector(left, right, sequencesVec[i], k ,rev,canon);
			fillSuffVector(left, right, sequencesVec[i], k ,rev,canon);
		}
	}
	//PARSE
	//seqsToRemoveInPref=seqsToRemoveInSuff={};
	parseVectorTip(left, right, sequencesVec, k, minSize);
	//TODO
	for (uint i(0); i < sequencesVec.size(); ++i){
		if (not sequencesVec[i].sequence.empty()){
			out << ">sequence_" + to_string(sequencesVec[i].index) << endl;
			out << sequencesVec[i].sequence << endl;
		}
	}
	auto endend=chrono::system_clock::now();auto waitedFor2=endend-startChrono;
	cout<<"whole process took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor2).count())<<" sec"<<endl;
	return 0;
}
