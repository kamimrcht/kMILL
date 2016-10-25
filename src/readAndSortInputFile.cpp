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

/* create buckets to store canonical sequences of the reads */
void createReadBuckets(uint nbBuckets, ifstream& readStructFile, vector <ofstream>& outFiles){
	string sequence, canonSequence;
	while (not readStructFile.eof()){
        getline(readStructFile, sequence);
		getline(readStructFile, sequence);
		canonSequence = getCanonical(sequence);
		if (not canonSequence.empty()){
			uint64_t key(transformStringToHash(canonSequence));
			outFiles[key % nbBuckets] << canonSequence << endl;
		}
	}
}


/* open buckets */
void openBuckets(vector<ofstream>& outFiles){
	for (uint nbFileOut(0); nbFileOut < outFiles.size() ; ++ nbFileOut){
		outFiles[nbFileOut].open("read_file_" + to_string(nbFileOut) + ".fa");
	}
}


/* get sequences from all buckets, sort them and remove duplicated reads */
void fillSortCleanBuckets(uint nbBuckets, vector <readStruct>& sequencesVec, uint threshold){
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
		if(threshold>0){
			cleanDuplicatesInreadStructs2(seqVecFile,threshold);
		}else{
			cleanDuplicatesInreadStructs(seqVecFile);
			}
		for (uint i(0); i < seqVecFile.size(); ++i){
			if (not seqVecFile[i].sequence.empty()){
				seqVecFile[i].index = index;
				++index;
				sequencesVec.push_back(seqVecFile[i]);
			}
		}
	}
}


/* remove bucket files */
void removeReadFiles(uint nbBuckets){
	for (uint i(0); i < nbBuckets; ++i){
		string s("read_file_" + to_string(i) + ".fa");
		remove(s.c_str());
	}
}


//test
void createReadBuckets2ndPass(uint nbBuckets, const vector <string>& sequences2ndPass, vector <ofstream>& outFiles){
	string sequence, canonSequence;
	for (uint i(0); i<sequences2ndPass.size(); ++i){
		canonSequence = getCanonical(sequences2ndPass[i]);
		if (not canonSequence.empty()){
			uint64_t key(transformStringToHash(canonSequence));
			outFiles[key % nbBuckets] << canonSequence << endl;
		}
	}
}
