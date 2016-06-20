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
#include "ograph.h"
#include <unordered_set>
#include <unordered_map>
#include <set>


using namespace std;



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./kMILL reads.fasta k" << endl;
		//~ srand (34567);
		//~ srand (time(NULL));
		//random reads
		//~ createinputlm(100*1000*1000,100);
		// random genome
		//~ randGenome(1000000);
		//genome from ref
		//~ perfectsReadsFromRef("simulGenome",100,1000*1000);
		//~ mutateReadsFromRef("../lambda_virus.fa", 100, 2*1000);
	} else {
		bool graph(false);
		if (argc == 4){
			string g(argv[3]);
			if (g == "-g"){
				graph = true;
			}
		}
		auto startChrono=chrono::system_clock::now();
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		ifstream readStructFile(fileName);
		uint nbBuckets(1);
		vector <ofstream> outFiles(nbBuckets);
		string titre("out_k"+to_string(k)+ "_" + getFileName(fileName) +".fa");
		ofstream out(titre);

		string sequence,sequence2;
		vector <readStruct> sequencesVec;
		openBuckets(outFiles);
		createReadBuckets(nbBuckets, readStructFile, outFiles);
		fillSortCleanBuckets(nbBuckets, sequencesVec);
		removeReadFiles(nbBuckets);
		setreadStructsIndex(sequencesVec);
		auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
		cout<<"Init took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
		unordered_set<string> seqsToRemoveInSuff;
		unordered_set<string> seqsToRemoveInPref;
		unordered_set<int> readsToRemovePref;
		unordered_set<int> readsToRemoveSuff;
		unordered_set<uint> colorNodePref;
		unordered_set<uint> colorNodeSuff;
		unordered_map<uint, uint> sizesNode;
		if (k>0) {
			do {
				auto startChrono=chrono::system_clock::now();
				vector <edge> right;  // vector of canonical suffixes
				vector <edge> left; //  vector of canonical prefixes
				auto startfor=chrono::system_clock::now();
				for (uint i(0); i<sequencesVec.size(); ++i){
					if (sequencesVec[i].sequence.size() >= k){
						fillPrefVector(left, right, sequencesVec[i], k, readsToRemovePref);
						fillSuffVector(left, right, sequencesVec[i], k, readsToRemoveSuff);
					}
				}
				auto endFor=chrono::system_clock::now();auto wFor=endFor-startfor;
				auto startparse=chrono::system_clock::now();
				if (graph){
					sequences2dot(sequencesVec, k, colorNodePref, colorNodeSuff, sizesNode);
					system("dot -Tpng out.dot > output.png");
					cin.get();
				}

				parseVector(left, right, sequencesVec, k, seqsToRemoveInSuff, seqsToRemoveInPref,  readsToRemovePref, readsToRemoveSuff);
				auto endparse=chrono::system_clock::now();auto wParse=endparse-startparse;
				auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
				cout<<"k: "<<k<<": left.size "<<left.size()<<" right.size "<<right.size()<<" Step took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec "<<endl;
				--k;				
			} while (k>2);
			if (graph){
				sequences2dot(sequencesVec, k, colorNodePref, colorNodeSuff, sizesNode);
				system("dot -Tpng out.dot > output.png");
			}
			for (uint i(0); i < sequencesVec.size(); ++i){
				if (not sequencesVec[i].sequence.empty()){
					out << ">sequence_" + to_string(sequencesVec[i].index) << endl;
					out << sequencesVec[i].sequence << endl;
				}
			}
			auto endend=chrono::system_clock::now();auto waitedFor2=endend-startChrono;
			cout<<"whole process took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor2).count())<<" sec"<<endl;
		}
	}
	return 0;
}
