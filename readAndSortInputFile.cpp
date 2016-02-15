#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "readAndSortInputFile.h"
#include "compaction.h"
#include "utils.h"


using namespace std;


hash<string> readHash;


uint64_t transformStringToHash(const string& read){
	return readHash(read);
}


void createReadBuckets(uint nbBuckets, ifstream& readStructFile, vector <ofstream>& outFiles){
	string sequence, canonSequence;
	while (not readStructFile.eof()){
		//~ cout<<1<<endl;
        getline(readStructFile, sequence);
		getline(readStructFile, sequence);
		canonSequence = getCanonical(sequence);
		if (not canonSequence.empty()){
			//~ cout<<sequence<<endl;
			uint64_t key(transformStringToHash(canonSequence));
			outFiles[key % nbBuckets] << canonSequence << endl;
		}
	}
}


void openBuckets(vector<ofstream>& outFiles){
	for (uint nbFileOut(0); nbFileOut < outFiles.size() ; ++ nbFileOut){
		outFiles[nbFileOut].open("read_file_" + to_string(nbFileOut) + ".fa");
	}
}


void fillSortCleanBuckets(uint nbBuckets, vector <readStruct>& sequencesVec){
	uint index(0);
	vector <readStruct> seqVecFile;
	string seq;
	for (uint nbFileOut(0); nbFileOut < nbBuckets ; ++ nbFileOut){
		seqVecFile={};
		ifstream readFile("read_file_" + to_string(nbFileOut) + ".fa");
		while (not readFile.eof()){
			getline(readFile, seq);
			if (not seq.empty()){
				seqVecFile.push_back({0, seq});
			}
		}
		sort(seqVecFile.begin(), seqVecFile.end(), compareRead());
		cleanDuplicatesInreadStructs(seqVecFile);
		for (uint i(0); i < seqVecFile.size(); ++i){
			if (not seqVecFile[i].sequence.empty()){
				seqVecFile[i].index = index;
				++index;
				sequencesVec.push_back(seqVecFile[i]);
			}
		}
	}
}


void removeReadFiles(uint nbBuckets){
	for (uint i(0); i < nbBuckets; ++i){
		string s("read_file_" + to_string(i) + ".fa");
		remove(s.c_str());
	}
}
