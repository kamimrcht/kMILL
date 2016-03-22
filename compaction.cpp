#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "utils.h"
#include "compaction.h"
#include "ograph.h"

#include <chrono>


using namespace std;


edge nPrefix(uint n, uint index, const string& sequence, bool canon){
	return {index, sequence.substr(0, n), canon};
}


edge nSuffix(uint n, uint index, const string& sequence, bool canon){
	return {index, sequence.substr(sequence.size()-n, n), canon};
}


vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff, uint k){
	vector<edge> vectResult;
	uint i(0);
	bool remove(false);
	while (i < vect.size()){
		if (i == 0){
			if (vect[i].sequence != vect[i+1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				remove = true;
			}
		} else if (i == vect.size()-1){
			if (vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				remove = true;
			}
		} else {
			if (vect[i].sequence != vect[i+1].sequence and vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				remove = true;
			}
		}
		if (remove == true){
			//~ cout << vect[i].sequence << endl;
			if (vect[i].canonical == true){
				seqsToRemoveInPref.insert(vect[i].sequence);
			} else {
				seqsToRemoveInSuff.insert(vect[i].sequence);
			}
		}
		++i;
		remove = false;
	}
	return vectResult;
}


vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff, uint k){
	vector<edge> vectResult;
	uint i(0);
	bool remove(false);
	while (i < vect.size()){
		if (i == 0){
			if (vect[i].sequence != vect[i+1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				//~ if (vect[i].sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
					//~ cout << "%%%%%%%%%%%%%%%%%%%%%%1" << endl;
				//~ }
				//~ cout << "cas 1"<< endl;
				remove = true;
			}
		} else if (i == vect.size()-1){
			if (vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				//~ if (vect[i].sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
					//~ cout << "%%%%%%%%%%%%%%%%%%%%%%2" << endl;
				//~ }
				//~ cout << "cas 2"<< endl;
				remove = true;
			}
		} else {
			if (vect[i].sequence != vect[i+1].sequence and vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			} else {
				//~ if (vect[i].sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
					//~ cout << "%%%%%%%%%%%%%%%%%%%%%%3" << endl;
				//~ }
				//~ cout << "cas 3"<< endl;
				remove = true;
			}
		}
		if (remove == true){
			//~ if (vect[i].sequence == "GTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
				//~ cout << "%%%%%%%%%%%%%%%%%%%%%%OUI" << endl;
			//~ }
			if (vect[i].canonical == true){
				//~ cout << "yepppppppppppppppppppp" << endl;
				//~ cout << vect[i].sequence<< " to rm in suff" << endl;
				seqsToRemoveInSuff.insert(vect[i].sequence);
			} else {
				seqsToRemoveInPref.insert(vect[i].sequence);
			}
		}
		++i;
		remove = false;
	}
	return vectResult;
}


void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove){
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		if (seqsToRemove.count(vectL[i].sequence)){
			if (vectL[i].canonical == true){  // prefix
				readsToRemove.insert(vectL[i].index);
			}
		}
		++i;
	}
	uint j(0);
	//~ cout << "%%%%%%%%%%%%%%%%%" << endl;
	//~ for (auto iter = seqsToRemove.begin(); iter != seqsToRemove.end(); ++iter){
		//~ cout << *iter << endl;
	//~ }
	while (j < vectR.size()){
		if (seqsToRemove.count(vectR[j].sequence)){
			if (vectR[j].canonical == false){  // prefix
				readsToRemove.insert(vectR[j].index);
			}
		}
		++j;
	}
	//~ cout << "------------------------" << endl;
	//~ for (auto iter = readsToRemove.begin(); iter != readsToRemove.end(); ++iter){
		//~ cout << *iter << endl;
	//~ }
	//~ cout << "%%%%%%%%%%%%%%%%%" << endl;
}


void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove){
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		if (seqsToRemove.count(vectL[i].sequence)){
			if (vectL[i].canonical == false){  // suffix
				readsToRemove.insert(vectL[i].index);
			}
		}
		++i;
	}
	uint j(0);
	//~ cout << "%%%%%%%%%%%%%%%%%" << endl;
	//~ for (auto iter = seqsToRemove.begin(); iter != seqsToRemove.end(); ++iter){
		//~ cout << *iter << endl;
	//~ }
	while (j < vectR.size()){
		if (seqsToRemove.count(vectR[j].sequence)){
			if (vectR[j].canonical == true){  // suffix
				readsToRemove.insert(vectR[j].index);
				//~ cout << vectR[j].index << endl;
			}
		}
		++j;
	}
	//~ cout << "------------------------" << endl;
	//~ for (auto iter = readsToRemove.begin(); iter != readsToRemove.end(); ++iter){
		//~ cout << *iter << endl;
	//~ }
	//~ cout << "%%%%%%%%%%%%%%%%%" << endl;
}



bool compareEdgeByString(const edge& seqL, const edge& seqR){
    return seqL.sequence < seqR.sequence;
}


//  compaction of two readStructs if they have a k-overlap
string compaction(const readStruct& seq1, const readStruct& seq2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
	
	edge beg1 = nPrefix(k, seq1.index, seq1.sequence, true);
	edge beg2 = nPrefix(k, seq2.index, seq2.sequence, true);
	edge end1 = nSuffix(k, seq1.index, seq1.sequence, true);
	edge end2 = nSuffix(k, seq2.index, seq2.sequence, true);
	edge rEnd1 = {end1.index, revComp(end1.sequence),true};
	edge rBeg2 = {beg2.index, revComp(beg2.sequence), true};
	string rSeq2 = revComp(seq2.sequence);
	if (end1.sequence == beg2.sequence){ //  overlap FF
		//~ cout << "éééééééééééééééééééééééééééééé " << seq1.index << " " << seq2.index << endl;
		string compacted(seq1.sequence + seq2.sequence.substr(k));
		string compactedToReturn(getCanonical(compacted));
		if (compacted != compactedToReturn){ // read's prefix will become its suffix and read's suffix will become its prefix
			if (readsToRemovePref.count(seq1.index)){
				//~ readsToRemovePref.erase(seq1.index);
				readsToRemoveSuff.insert(seq1.index);
			}
			if (readsToRemoveSuff.count(seq2.index)){
				//~ readsToRemoveSuff.erase(seq2.index);
				readsToRemovePref.insert(seq1.index);
			}
		} else { // register read's new suffix (the sequence 2's suffix)
			readsToRemoveSuff.insert(seq1.index);
		}	
		return compactedToReturn;
	} else if (end2.sequence == beg1.sequence) { //  overlap RR
		string compacted(seq2.sequence + seq1.sequence.substr(k));
		string compactedToReturn(getCanonical(compacted));
		if (compacted != compactedToReturn){ 
			if (readsToRemoveSuff.count(seq1.index)){
				//~ readsToRemoveSuff.erase(seq1.index);
				readsToRemovePref.insert(seq1.index);
			}
			if (readsToRemovePref.count(seq2.index)){
				//~ readsToRemovePref.erase(seq2.index);
				readsToRemoveSuff.insert(seq1.index);
			}
		} else { // register read's new prefix
			readsToRemovePref.erase(seq2.index);
		}
		return compactedToReturn;
	} else if (rEnd1.sequence == end2.sequence) { //  overlap FR
		string compacted(seq1.sequence + rSeq2.substr(k));
		string compactedToReturn(getCanonical(compacted));
		if (compacted == compactedToReturn){ // register read's new suffix which is read 2's prefix
			if (readsToRemovePref.count(seq2.index)){
				//~ readsToRemovePref.erase(seq2.index);
				readsToRemoveSuff.insert(seq1.index);
			}
		} else { // new suffix is the read 1's prefix, new prefix is read 2's prefix
			if (readsToRemovePref.count(seq1.index)){
				//~ readsToRemovePref.erase(seq1.index);
				readsToRemoveSuff.insert(seq1.index);
			}
			if (readsToRemovePref.count(seq2.index)){
				//~ readsToRemovePref.erase(seq2.index);
				readsToRemovePref.insert(seq1.index);
			}
		}
		//~ return seq1.sequence + rSeq2.substr(k);
		return compactedToReturn;
	} else if (beg1.sequence == rBeg2.sequence){ //  overlap RF
		//~ cout << "lolol" << endl;
		string compacted(rSeq2 + seq1.sequence.substr(k));
		string compactedToReturn(getCanonical(compacted));
		if (compacted == compactedToReturn){ // register read's new prefix which is read 2's suffix
			if (readsToRemoveSuff.count(seq2.index)){
				readsToRemoveSuff.erase(seq2.index);
				readsToRemovePref.insert(seq1.index);
			}
		} else { // new read's suffix is read 2's suffix, new read's prefix is read 1's suffix
			if (readsToRemoveSuff.count(seq1.index)){
				readsToRemoveSuff.erase(seq1.index);
				readsToRemovePref.insert(seq1.index);
			}
			if (readsToRemoveSuff.count(seq2.index)){
				readsToRemoveSuff.erase(seq2.index);
				readsToRemoveSuff.insert(seq1.index);
			}
		}
		return compactedToReturn;
		//~ return rSeq2 + seq1.sequence.substr(k);
	} else {
		cout<<"fail..."<<endl;
		return "";
	}
}


//  If two readStructs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate.
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
	if (not vec[indexreadStruct1].sequence.empty()){
		if (not vec[indexreadStruct2].sequence.empty()){
			string c = compaction(vec[indexreadStruct1], vec[indexreadStruct2], k, readsToRemovePref, readsToRemoveSuff);
			/* debug */
			bool present(isSubSequenceInSequence("AGCGCGGTGGTCCCACCTGACCCCATGCCGAACTCAGAAGTGAAACGCCGTAGCGCCGATGGTAGTGTGGGGTCTCCCCATGCGAGAGTAGGGAACTGCC", c));
			if (present == true){
				cout << "compac " << c << " "<<k << " " << vec[indexreadStruct1].index << " " <<vec[indexreadStruct1].sequence << " " << vec[indexreadStruct2].index << " "<< vec[indexreadStruct2].sequence<<endl;
				cin.get();
				//~ if (readsToRemovePref.count(indexreadStruct1)){
					//~ cout << "pref read 1" << endl;
				//~ }
				//~ if (readsToRemoveSuff.count(indexreadStruct1)){
					//~ cout << "suff read 1" << endl;
				//~ }
				//~ if(readsToRemovePref.count(indexreadStruct2)){
					//~ cout << "pref read 2" << endl;
				//~ }
				//~ if(readsToRemoveSuff.count(indexreadStruct2)){
					//~ cout << "suff read 2" << endl;
				//~ }
			}
			/* end db */
			if (not c.empty()){
				//~ cout << "COMPACTIOn "  << indexreadStruct1 << " " <<vec[indexreadStruct1].sequence << " " << indexreadStruct2 << " " << vec[indexreadStruct2].sequence<<endl;
				vec[indexreadStruct1] = {vec[indexreadStruct1].index, c};
				vec[indexreadStruct2].index = vec[indexreadStruct1].index;
				vec[indexreadStruct2].sequence = "";
			}
		} else {
			compactInVector(vec, indexreadStruct1, vec[indexreadStruct2].index, k,  readsToRemovePref, readsToRemoveSuff); //  each time a sequence is empty, the index leads to the sequence it's been compacted in-> recursive call until we find the sequence
		}
	} else {
		if (not vec[indexreadStruct2].sequence.empty()){
			compactInVector(vec, vec[indexreadStruct1].index, indexreadStruct2, k, readsToRemovePref, readsToRemoveSuff);
		} else {
			compactInVector(vec, vec[indexreadStruct1].index, vec[indexreadStruct2].index, k, readsToRemovePref, readsToRemoveSuff);
		}
	}
}


//  checks from the suffixes and prefixes of pairs of readStructs of a vector if they can be compacted
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemoveInSuff, unordered_set<string>& seqsToRemoveInPref, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
	auto startrmappend=chrono::system_clock::now();
	uint compac(0);
	sort(left.begin(), left.end(), compareEdge());
	sort(right.begin(), right.end(), compareEdge());
	vector<edge> leftRemoved;
	vector<edge> rightRemoved;
	vector<edge> leftSingles;
	vector<edge> rightSingles;
	if (left.size()>1){
		leftSingles = removeNotSinglesInLeft(left, seqsToRemoveInPref, seqsToRemoveInSuff, k);
	} else {
		leftSingles = left;
	}
	if (right.size()>1){
		rightSingles = removeNotSinglesInRight(right, seqsToRemoveInPref, seqsToRemoveInSuff, k);
	} else {
		rightSingles = right;
	}
	//~ for (uint i(0); i < leftSingles.size(); ++i){
		//~ cout <<"L " <<leftSingles[i].sequence <<leftSingles[i].index << endl;
	//~ }
	//~ for (uint i(0); i < rightSingles.size(); ++i){
		//~ cout << "R "<< rightSingles[i].sequence << rightSingles[i].index << endl;
	//~ }
	appendListReadsToRemovePref(leftSingles, rightSingles, seqsToRemoveInSuff, readsToRemovePref);
	appendListReadsToRemoveSuff(leftSingles, rightSingles, seqsToRemoveInPref, readsToRemoveSuff);
	int mmm = 4461;
	if (readsToRemoveSuff.count(mmm)==1){
		cout << "IN SUFF" << endl;
	}
	if (readsToRemovePref.count(mmm)==1){
		cout << "IN PREF" << endl;
	}
	/* db */
	appendListReadsToRemovePref(leftSingles, rightSingles, seqsToRemoveInSuff, readsToRemoveSuff);
	appendListReadsToRemoveSuff(leftSingles, rightSingles, seqsToRemoveInPref, readsToRemovePref);
	auto endrmappend=chrono::system_clock::now();auto wrmappend=endrmappend-startrmappend;
	cout<<"parse 1 took : "<<(chrono::duration_cast<chrono::seconds>(wrmappend).count())<<" sec"<<endl;
	/* end db*/
	auto startleft=chrono::system_clock::now();
	uint indexL(0),indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			if (leftSingles[indexL].index != rightSingles[indexR].index){
				compactInVector(readStructsVec, leftSingles[indexL].index, rightSingles[indexR].index, k,  readsToRemovePref, readsToRemoveSuff);
				//~ cout << "compaction " << k << endl;
				++compac;
			}
			++indexL;
			++indexR;
		} else {
			if (leftSingles[indexL].sequence < rightSingles[indexR].sequence){
				++indexL;
			} else {
				++indexR;
			}
		}
	}
	auto endleft=chrono::system_clock::now();auto wleft=endleft-startleft;
	cout<<"parse 1 took : "<<(chrono::duration_cast<chrono::seconds>(wleft).count())<<" sec"<<endl;
}


// fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemovePref){
	if (not readsToRemovePref.count(seq.index)){
		edge prefix = nPrefix(k, seq.index, seq.sequence, true);
		string canonPrefix = getCanonical(prefix.sequence);
		string s(getCanonical("GTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"));

		if (prefix.sequence == canonPrefix){
			//~ if (prefix.sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
			//~ if (prefix.sequence ==	s){
				//~ cout << "in L" << endl;
			//~ }
			vecLeft.push_back({prefix.index, canonPrefix, true});
		} else {
			//~ if (prefix.sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"){
			//~ if (prefix.sequence ==	s){
				//~ cout << "in R" << endl;
			//~ }
			vecRight.push_back({prefix.index, canonPrefix, false});
		}
	}
}


// fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemoveSuff){
	if (not readsToRemoveSuff.count(seq.index)){
		edge suffix = nSuffix(k, seq.index, seq.sequence, true);
		
		string canonSuffix = getCanonical(suffix.sequence);
		string s(getCanonical("GTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCA"));
		if (suffix.sequence == canonSuffix){
			
			//~ if (suffix.sequence ==	s){
				//~ cout << "in R" << endl;
			//~ }
			vecRight.push_back({suffix.index, canonSuffix, true});
		} else {
			//~ if (suffix.sequence == "CGGTGGCGGCGGCACAGTTGGTGTGGCGGCCTCAGTCCGGAACAATTTGAAAACAAGAACCTCGCTTAGGCCTGTGTCCATATTACGTGGGTAGGATCAA"){
			//~ if (suffix.sequence ==	s){
				//~ cout << "in L" << endl;
			//~ }
			vecLeft.push_back({suffix.index, canonSuffix, false});
		}
	}
}


void cleanDuplicatesInreadStructs(vector <readStruct>& vec){
	uint i(0);
	string previousSeq("");
	while(i<vec.size()){
		string temp = vec[i].sequence;
		if (temp == previousSeq){
			vec[i].sequence = "";
		}else{
			previousSeq = temp;
		}
		++i;
	}
}


void setreadStructsIndex(vector <readStruct>& vec){
	for (uint i(0); i<vec.size(); ++i){
		if (not vec[i].sequence.empty()){
			vec[i].index = i;
		}
	}
}


void initVectofreadStructs(vector <readStruct>& vec, const string& sequence){
	if (not sequence.empty()){
		vec.push_back({0, sequence});
	}
}
