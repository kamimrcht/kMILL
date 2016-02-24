#include <fstream>
#include <cstring>
#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include "compaction.h"
#include "readAndSortInputFile.h"
#include "utils.h"

using namespace std;



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./kMILL reads.fasta k" << endl;
	} else {
		// srand (34567);
		srand (time(NULL));
		//random genome
		// createinputlm(2*1000*1000,100);
		//genome from ref
		// prefectsReadsFromRef("ecoliref.fa",100,1*1000*1000);
		auto startChrono=chrono::system_clock::now();
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		ifstream readStructFile(fileName);
		uint nbBuckets(1);
		vector <ofstream> outFiles(nbBuckets);
		ofstream out("out_k"+to_string(k)+ "_" + fileName +".fa");
		string sequence,sequence2;
		vector <readStruct> sequencesVec;
		openBuckets(outFiles);
		createReadBuckets(nbBuckets, readStructFile, outFiles);
		fillSortCleanBuckets(nbBuckets, sequencesVec);
		removeReadFiles(nbBuckets);
		setreadStructsIndex(sequencesVec);
		auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
		cout<<"Init took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
		do {
			auto startChrono=chrono::system_clock::now();

			vector <edge> right;  // vector of canonical suffixes
			vector <edge> left; //  vector of canonical prefixes
			for (uint i(0); i<sequencesVec.size(); ++i){
				if (sequencesVec[i].sequence.size() > k){
					fillPrefVector(left, right, sequencesVec[i], k);
					fillSuffVector(left, right, sequencesVec[i], k);
				}
			}

			parseVector(left, right, sequencesVec, k);
			//~ parseVector(left, right, sequencesVec, k, outF);
			--k;
			auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
			cout<<"k: "<<k<<": left.size "<<left.size()<<" right.size "<<right.size()<<" Step took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec "<<endl;
		} while (k>2);
		uint header(1);
		for (uint i(0); i<sequencesVec.size(); ++i){
			if (not sequencesVec[i].sequence.empty()){
				out<< ">seq_" + to_string(header) << endl;
				out<<sequencesVec[i].sequence << endl;
				++header;
			}
		}
		auto endend=chrono::system_clock::now();auto waitedFor2=endend-startChrono;
		cout<<"whole process took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor2).count())<<" sec"<<endl;
	}
	return 0;
}
