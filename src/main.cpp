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


using namespace std;



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./kMILL reads.fasta k" << endl;
		//~ srand (34567);
		//~ srand (time(NULL));
		//random reads
		//~ createinputlm(100*1000*1000,100);
		// random genome
		randGenome(1000);
		//genome from ref
		perfectsReadsFromRef("simulGenome",100,100);
		//~ perfectsReadsFromRef("../ecoliref.fa",200,1*1000*1000);
		//~ perfectsReadsFromRef("../lambda_virus.fa",100,5*100*1000);
		//~ perfectsReadsFromRef("simulGenome",200,1*1000*1000);
	} else {
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
		//~ uint size(sequencesVec.size());
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
				/* debug */
				if (k == 94){
					for (uint i(0); i < sequencesVec.size(); ++i){
						if (sequencesVec[i].index == 52){
							cout << k << " " << sequencesVec[i].index << " "<< sequencesVec[i].sequence << " "<< sequencesVec[i].sequence.size() << endl;
						}
					}
				}
				if (k == 93){
					for (uint i(0); i < sequencesVec.size(); ++i){
						if (sequencesVec[i].index == 52){
							cout << k << " " <<sequencesVec[i].index << " "<< sequencesVec[i].sequence <<" "<< sequencesVec[i].sequence.size() <<endl;
						}
					}
				}
				sequences2dot(sequencesVec, k, colorNodePref, colorNodeSuff, sizesNode);
				system("dot -Tpng out.dot > output.png");
				cin.get();
				/* debug */
				parseVector(left, right, sequencesVec, k, seqsToRemoveInSuff, seqsToRemoveInPref,  readsToRemovePref, readsToRemoveSuff);
				auto endparse=chrono::system_clock::now();auto wParse=endparse-startparse;
				auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
				cout<<"k: "<<k<<": left.size "<<left.size()<<" right.size "<<right.size()<<" Step took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec "<<endl;
				/* debug */
				//~ ifstream readStructFile2(fileName);
				//~ unordered_set <string> initSet;
				//~ unordered_map <string, uint> finalSet;
				//~ uint inserted(0);
				//~ while (not readStructFile2.eof()){
					//~ string seq;
					//~ getline(readStructFile2, seq);
					//~ getline(readStructFile2, seq);
					//~ string seqC(getCanonical(seq));
					//~ if (not seqC.empty()){
						//~ initSet.insert(seqC);
						//~ ++inserted;
					//~ }
				//~ }
				//~ uint rfound(0);
				//~ for (uint i(0); i < sequencesVec.size(); ++i){
					//~ if (not sequencesVec[i].sequence.empty()){
						//~ uint w(0);
						//~ do{
							//~ readStruct read;
							//~ read.index = sequencesVec[i].index;
							//~ read.sequence = getCanonical(sequencesVec[i].sequence.substr(w,100));
							//~ ++w;
							//~ if (not read.sequence.empty()){
								//~ auto readFound = initSet.find(read.sequence);
								//~ if (readFound == initSet.end() ){
								//~ } else {
									//~ ++rfound;
									//~ if (finalSet.count(read.sequence)){
										//~ finalSet[read.sequence] += 1;
											//~ cout << read.index << " REPETE: " <<read.sequence <<  " at "<< k <<  endl;
									//~ } else {
										//~ finalSet.insert({read.sequence,1});
									//~ }
								//~ }
							//~ }

						//~ } while(w < sequencesVec[i].sequence.size());
					//~ }
				//~ }
				//~ if ( initSet.size() == finalSet.size()){
					//~ cout << "*** debug *** all reads in unitigs: ok" << endl;
					//~ cout <<  "sizes of sets: " << initSet.size() << " " << finalSet.size() << endl;
					//~ cout <<  "inserted elements: " <<inserted << " " << rfound << endl;
				//~ } else {
					//~ cout << "size real reads set: " <<  size << endl;
					//~ cout <<  "sizes of sets: " << initSet.size() << " " << finalSet.size() << endl;
					//~ cout <<  "inserted elements: " <<inserted << " " << rfound << endl;
				//~ }
				/* end debug*/

				

				--k;				
			} while (k>2);
			sequences2dot(sequencesVec, k, colorNodePref, colorNodeSuff, sizesNode);
			system("dot -Tpng out.dot > output.png");
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
