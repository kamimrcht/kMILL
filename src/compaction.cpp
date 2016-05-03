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

/* remove overlaps that are duplicated in vector right, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
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


/* remove overlaps that are duplicated in vector left, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff, uint k){
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
			if (vect[i].canonical == true){
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


/* from sequences spotted, get reads' index for reads we do not want the prefix to be present in the next pass */ 
void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove){
	//~ cout << "rmv prf" << endl;
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		//~ cout << "size pour mon mari "<< vectL[i].sequence<<endl;
		//~ cout << "size du sex " << seqsToRemove.size() << endl;
		if (seqsToRemove.unordered_set::count(vectL[i].sequence)){
			//~ cout << "1" << endl;
			if (vectL[i].canonical == true){  // prefix
				//~ cout << "11" << endl;
				readsToRemove.insert(vectL[i].index);
			}
		}
		++i;
	}
	uint j(0);
	while (j < vectR.size()){
		if (seqsToRemove.unordered_set::count(vectR[j].sequence)){
			//~ cout << "2" << endl;
			if (vectR[j].canonical == false){  // prefix
				readsToRemove.insert(vectR[j].index);
			}
		}
		++j;
	}
}


/* from sequences spotted, get reads' index for reads we do not want the suffix to be present in the next pass */ 
void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, unordered_set<int>& readsToRemove){
	//~ cout << "rmv sff" << endl;
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		//~ cout << "size pour mon mari "<< vectL[i].sequence<<endl;
		//~ cout << "size du sex " << seqsToRemove.size() << endl;
		if (seqsToRemove.unordered_set::count(vectL[i].sequence)){
			//~ cout << "3" << endl;
			if (vectL[i].canonical == false){  // suffix
				readsToRemove.insert(vectL[i].index);
			}
		}
		++i;
	}
	uint j(0);
	while (j < vectR.size()){
		if (seqsToRemove.unordered_set::count(vectR[j].sequence)){
			//~ cout << "4" << endl;
			if (vectR[j].canonical == true){  // suffix
				readsToRemove.insert(vectR[j].index);
			}
		}
		++j;
	}
}



bool compareEdgeByString(const edge& seqL, const edge& seqR){
    return seqL.sequence < seqR.sequence;
}


/*  compaction of two unitigs if they have a k-overlap */
string compaction(const readStruct& seq1, const readStruct& seq2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
	
	edge beg1 = nPrefix(k, seq1.index, seq1.sequence, true);
	edge beg2 = nPrefix(k, seq2.index, seq2.sequence, true);
	edge end1 = nSuffix(k, seq1.index, seq1.sequence, true);
	edge end2 = nSuffix(k, seq2.index, seq2.sequence, true);
	edge rEnd1 = {end1.index, revComp(end1.sequence), true};
	edge rBeg2 = {beg2.index, revComp(beg2.sequence), true};
	string rSeq2 = revComp(seq2.sequence);
	if (end1.sequence == beg2.sequence){ //  overlap FF
		string compacted(seq1.sequence + seq2.sequence.substr(k));
		bool b(isCanonical(compacted));
		if (not b){ // read's prefix will become its suffix and read's suffix will become its prefix
			string compactedToReturn(getCanonical(compacted));
			if (readsToRemovePref.unordered_set::count(seq1.index)){
				//~ cout <<  "AAAAAAAAAAAA" << endl;
				readsToRemoveSuff.insert(seq1.index);
			}
			if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				//~ cout << "BBBBBBBBBBBBBB" << endl;
				readsToRemovePref.insert(seq1.index);
			}
			return compactedToReturn;
		} else { // register read's new suffix (the sequence 2's suffix)
			//~ cout << "CCCCCCCCCCCC" << endl;
			if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				readsToRemoveSuff.insert(seq1.index);
			}
			return compacted;
		}
	} else if (end2.sequence == beg1.sequence) { //  overlap RR
		string compacted(seq2.sequence + seq1.sequence.substr(k));
		bool b(isCanonical(compacted));
		if (not b){
			string compactedToReturn(getCanonical(compacted));
			if (readsToRemoveSuff.unordered_set::count(seq1.index)){
				readsToRemovePref.insert(seq1.index);
			}
			if (readsToRemovePref.unordered_set::count(seq2.index)){
				readsToRemoveSuff.insert(seq1.index);
			}
			return compactedToReturn;
		} else { // register read's new prefix
			//~ readsToRemovePref.erase(seq2.index);
			return compacted;
		}
	} else if (rEnd1.sequence == end2.sequence) { //  overlap FR
		string compacted(seq1.sequence + rSeq2.substr(k));
		bool b(isCanonical(compacted));
		if (not b){ // new suffix is the read 1's prefix, new prefix is read 2's prefix
			string compactedToReturn(getCanonical(compacted));
			if (readsToRemovePref.unordered_set::count(seq1.index)){
				readsToRemoveSuff.insert(seq1.index);
			}
			if (readsToRemovePref.unordered_set::count(seq2.index)){
				readsToRemovePref.insert(seq1.index);
			}
			return compactedToReturn;
		} else { // register read's new suffix which is read 2's prefix
			if (readsToRemovePref.unordered_set::count(seq2.index)){
				readsToRemoveSuff.insert(seq1.index);
			}
			return compacted;
		}
	} else if (beg1.sequence == rBeg2.sequence){ //  overlap RF
		string compacted(rSeq2 + seq1.sequence.substr(k));
		/* tests */
		//~ string compactedToReturn(getCanonical(compacted));
		bool b(isCanonical(compacted));
		if (not b){ // new read's suffix is read 2's suffix, new read's prefix is read 1's suffix
			string compactedToReturn(getCanonical(compacted));
			if (readsToRemoveSuff.unordered_set::count(seq1.index)){
				readsToRemovePref.insert(seq1.index);
			}
			if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				readsToRemoveSuff.insert(seq1.index);
			}
			return compactedToReturn;
		} else { // register read's new prefix which is read 2's suffix
			if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				readsToRemovePref.insert(seq1.index);
			}
			return compacted;
		}
	} else {
		cout<<"fail..."<<endl;
		return "";
	}
}


/*  Recursive func to compact unitigs (if they are sequences or indexes)). If two unitigs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate. */
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
	if (not vec[indexreadStruct1].sequence.empty()){
		if (not vec[indexreadStruct2].sequence.empty()){
			string c = compaction(vec[indexreadStruct1], vec[indexreadStruct2], k, readsToRemovePref, readsToRemoveSuff);
			if (not c.empty()){
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


/*  checks from the suffixes and prefixes of pairs of unitigs if they can be compacted. Appends overlaps that should be removed in the next pass. */
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemoveInSuff, unordered_set<string>& seqsToRemoveInPref, unordered_set<int>& readsToRemovePref, unordered_set<int>& readsToRemoveSuff){
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
	//~ cout << "SIZE " << seqsToRemoveInPref.size() << endl;
	//~ cout << "SIZE " <<  seqsToRemoveInSuff.size() << endl;
	appendListReadsToRemovePref(leftSingles, rightSingles, seqsToRemoveInPref, readsToRemovePref);
	appendListReadsToRemoveSuff(leftSingles, rightSingles, seqsToRemoveInSuff, readsToRemoveSuff);
	//~ cout << "SIZE " << readsToRemovePref.size() << endl;
	//~ cout << "SIZE " <<  readsToRemoveSuff.size() << endl;
	uint indexL(0),indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			if (leftSingles[indexL].index != rightSingles[indexR].index){
				compactInVector(readStructsVec, leftSingles[indexL].index, rightSingles[indexR].index, k,  readsToRemovePref, readsToRemoveSuff);
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
	//~ cout << "SIZE " << readsToRemovePref.size() << endl;
	//~ cout << "SIZE " <<  readsToRemoveSuff.size() << endl;
}


/* fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs */
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemovePref){
	if (not readsToRemovePref.unordered_set::count(seq.index)){
		edge prefix = nPrefix(k, seq.index, seq.sequence, true);
		string canonPrefix = getCanonical(prefix.sequence);
		if (prefix.sequence == canonPrefix){
			vecLeft.push_back({prefix.index, canonPrefix, true});
		} else {
			vecRight.push_back({prefix.index, canonPrefix, false});
		}
	}
}


/* fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs */
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k, const unordered_set<int>& readsToRemoveSuff){
	if (not readsToRemoveSuff.unordered_set::count(seq.index)){
		edge suffix = nSuffix(k, seq.index, seq.sequence, true);
		string canonSuffix = getCanonical(suffix.sequence);
		if (suffix.sequence == canonSuffix){
			vecRight.push_back({suffix.index, canonSuffix, true});
		} else {
			vecLeft.push_back({suffix.index, canonSuffix, false});
		}
	}
}


/* remove duplicates in reads*/
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


/* give a proper index to reads */
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
