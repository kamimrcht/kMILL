#include <fstream>
#include <cstring>
#include <iostream>
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
	} else {
		//~ srand (34567);
		srand (time(NULL));
		//random genome
		//~ createinputlm(2*1000*1000,100);
		//genome from ref
		//~ perfectsReadsFromRef("ecoliref.fa",100,1*1000*1000);
		//~ perfectsReadsFromRef("lambda_virus.fa",100,1*1000*5);
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
		uint size(sequencesVec.size());
		unordered_set<string> edgesToRemove;
		if (k>0) {
			do {
				bool diff = false;
				graph3 graphCompactor(k,0,1,sequencesVec.size());
				auto startChrono=chrono::system_clock::now();
				vector <edge> right;  // vector of canonical suffixes
				vector <edge> left; //  vector of canonical prefixes
				for (uint i(0); i<sequencesVec.size(); ++i){
					if (sequencesVec[i].sequence.size() >= k){
						graphCompactor.addtuple(make_tuple(sequencesVec[i].sequence,0,0));
						fillPrefVector(left, right, sequencesVec[i], k);
						fillSuffVector(left, right, sequencesVec[i], k);
					}
				}
				graphCompactor.debruijn();
				parseVector(left, right, sequencesVec, k, edgesToRemove);
				unordered_set <string> malfoySet;
				unordered_set <string> camilleSet;
				uint sizeM(0);

				for (uint i(0); i < graphCompactor.size(); ++i){
					string seq = getCanonical(graphCompactor.unitigs[i]);
					if (not seq.empty()){
						malfoySet.insert(seq);
						sizeM += seq.size();
					}
				}

				uint sizeC(0);
				for (uint i(0); i<sequencesVec.size(); ++i){
					string seq = getCanonical(sequencesVec[i].sequence);
					sizeC += seq.size();
					if (not (seq.empty()) and  seq.size()>= k){
						camilleSet.insert(seq);
						if((malfoySet.count(seq)==0) ){
							cout << "*** c and not m ***" << seq << " "<<  seq.size() <<endl;
							diff = true;
						}
					}
				}
				for (uint i(0); i < graphCompactor.size(); ++i){
					string seq = getCanonical(graphCompactor.unitigs[i]);
					if (not seq.empty()){
						if(camilleSet.count(seq)==0){
							diff=true;
						}
					}
				}
				cout <<  "size gphe " << graphCompactor.size() << endl;
				cout <<"m:"<< sizeM << " c: " << sizeC << endl;
				auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
				cout<<"k: "<<k<<": left.size "<<left.size()<<" right.size "<<right.size()<<" Step took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec "<<endl;
				--k;
			} while (k>3);




			/* debug */
			ifstream readStructFile2(fileName);
			unordered_set <string> initSet;
			unordered_map <string, uint> finalSet;
			uint inserted(0);
			while (not readStructFile2.eof()){
				string seq;
				getline(readStructFile2, seq);
				getline(readStructFile2, seq);
				string seqC(getCanonical(seq));
				if (not seqC.empty()){
					initSet.insert(seqC);
					++inserted;
				}
			}
			uint rfound(0);
			for (uint i(0); i < sequencesVec.size(); ++i){
				if (not sequencesVec[i].sequence.empty()){
					uint w(0);
					do{
						readStruct read;
						if (w +100 < sequencesVec[i].sequence.size()) {
							read.index = sequencesVec[i].index;
							read.sequence = getCanonical(sequencesVec[i].sequence.substr(w,100));
						} else {
							read.index = sequencesVec[i].index;
							read.sequence = getCanonical(sequencesVec[i].sequence.substr(w));
						}

						++w;
						if (not read.sequence.empty()){
							auto readFound = initSet.find(read.sequence);
							if (readFound == initSet.end() ){
							} else {
								++rfound;
								if (finalSet.count(read.sequence)){
									finalSet[read.sequence] += 1;
										cout << read.index << " REPETE: " <<read.sequence << endl;
								} else {
									finalSet.insert({read.sequence,1});
								}
							}
						}

					} while(w < sequencesVec[i].sequence.size());
				}
			}
			//~ AGCACCGGTGATCATGTTTTTAACATAGTCGGCGTGCCCCGGGCAGTCTACGTGTGCGTAGTGACGGGTCGGGGTGTCGTATTCAACGTGAGAAGTGTTG
			//~ CAACACTTCTCACGTTGAATACGACACCCCGACCCGTCACTACGCACACGTAGACTGCCCGGGGCACGCCGACTATGTTAAAAACATGATCACCGGTGCT
			/* end debug*/
			uint header(1);
			vector <readStruct> sequencesVec2ndPass;
			sequencesVec2ndPass = sequencesVec;
			//~ uint gg(0);
			for (uint i(0); i<sequencesVec.size(); ++i){
				if (not sequencesVec[i].sequence.empty()){
					out<< ">seq_" + to_string(sequencesVec[i].index) << endl;
					out<<sequencesVec[i].sequence << endl;
					++header;
				}
			}
			//~ // test : 2nd pass
			//~ k = stoi(argv[2]);
			//~ ofstream out2nd("out_k"+to_string(k)+ "_2pass_" + fileName +".fa");
			//~ do {
				//~ auto startChrono=chrono::system_clock::now();
				//~ vector <edge> right;  // vector of canonical suffixes
				//~ vector <edge> left; //  vector of canonical prefixes
				//~ for (uint i(0); i<sequencesVec2ndPass.size(); ++i){
					//~ if (sequencesVec2ndPass[i].sequence.size() > k){
						//~ fillPrefVector(left, right, sequencesVec2ndPass[i], k);
						//~ fillSuffVector(left, right, sequencesVec2ndPass[i], k);
					//~ }
				//~ }

				//~ parseVector(left, right, sequencesVec2ndPass, k);
				//~ cout<<"k: "<<k<<": left.size "<<left.size()<<" right.size "<<right.size()<<" Step took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec "<<endl;
				//~ --k;
				//~ auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
			//~ } while (k>2);
			//~ uint header2(1);
			//~ for (uint i(0); i<sequencesVec2ndPass.size(); ++i){
				//~ if (not sequencesVec2ndPass[i].sequence.empty()){
					//~ out2nd<< ">seq_" + to_string(header2) << endl;
					//~ out2nd<<sequencesVec2ndPass[i].sequence << endl;
					//~ ++header;

				//~ }
			//~ }

			auto endend=chrono::system_clock::now();auto waitedFor2=endend-startChrono;
			cout<<"whole process took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor2).count())<<" sec"<<endl;
			/* debug */
			if ( initSet.size() == finalSet.size()){
				cout << "*** debug *** all reads in unitigs: ok" << endl;
				cout <<  "sizes of sets: " << initSet.size() << " " << finalSet.size() << endl;
				cout <<  "inserted elements: " <<inserted << " " << rfound << endl;
			} else {
				cout << "size real reads set: " <<  size << endl;
				cout <<  "sizes of sets: " << initSet.size() << " " << finalSet.size() << endl;
				cout <<  "inserted elements: " <<inserted << " " << rfound << endl;
			}
		}
	}
	return 0;
}
