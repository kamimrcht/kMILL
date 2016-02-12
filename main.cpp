#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "compaction.h"
#include "readAndSortInputFile.h"

using namespace std;



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./kMILL reads.fasta k" << endl;
	} else {
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		ifstream readStructFile(fileName);
		uint nbBuckets(10);
		vector <ofstream> outFiles(nbBuckets);
		ofstream out("out.fa");
		string sequence,sequence2;
		vector <readStruct> sequencesVec;
		openBuckets(outFiles);
		createsReadBuckets(nbBuckets, readStructFile, outFiles);
		fillSortCleanBuckets(nbBuckets, sequencesVec);
		removeReadFiles(nbBuckets);
		setreadStructsIndex(sequencesVec);
		do {
            cout<<k<<" Please be patient..."<<endl;
			vector <edge> right;  // vector of canonical suffixes
			vector <edge> left; //  vector of canonical prefixes
			for (uint i(0); i<sequencesVec.size(); ++i){
				if (sequencesVec[i].sequence.size() > k){
					fillPrefVector(left, right, sequencesVec[i], k);
					fillSuffVector(left, right, sequencesVec[i], k);
				}
			}
			parseVector(left, right, sequencesVec, k);
			--k;
		} while (k>2);
		for (uint i(0); i<sequencesVec.size(); ++i){
			if (not sequencesVec[i].sequence.empty()){
				out<<sequencesVec[i].sequence << endl;
			}
		}
	}
	return 0;
}
