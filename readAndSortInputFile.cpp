#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "readAndSortInputFile.h"
#include "compaction.h"

using namespace std;


uint64_t transformStringToHash(string read){
	hash<string> readHash;
	return readHash(read);
}


void createReadBuckets(uint nbBuckets, ifstream& readStructFile, vector <ofstream>& outFiles){
	while (not readStructFile.eof()){
		//~ cout<<1<<endl;
		string sequence;
        getline(readStructFile, sequence);
		getline(readStructFile, sequence);
		if (not sequence.empty()){
			//~ cout<<sequence<<endl;
			uint64_t key(transformStringToHash(sequence));
			int numFile(key % nbBuckets);
			outFiles[numFile] << sequence << endl;
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
	for (uint nbFileOut(0); nbFileOut < nbBuckets ; ++ nbFileOut){
		vector <readStruct> seqVecFile;
		string seq;
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
