#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "utils.h"
#include <chrono>



using namespace std;
uint nBucketsOverlap(100);



hash<std::string> strHash;



edge nPrefix(uint n, uint index, const string& sequence, bool canon){
	return {index, sequence.substr(0, n), canon};
}



edge nSuffix(uint n, uint index, const string& sequence, bool canon){
	return {index, sequence.substr(sequence.size()-n, n), canon};
}



/* remove overlaps that are duplicated in vector right, and spot sequences that should be remove in next pass because they reached their best overlap with a duplicated sequence */
vector<edge> removeNotSinglesInRight(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff){
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
vector<edge> removeNotSinglesInLeft(const vector<edge>& vect, unordered_set<string>& seqsToRemoveInPref, unordered_set<string>& seqsToRemoveInSuff){
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
void appendListReadsToRemovePref(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, vector<readStruct>& readStructsVec){
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		if (seqsToRemove.unordered_set::count(vectL[i].sequence)){
			if (vectL[i].canonical){  // prefix
				//~ readsToRemove.insert(vectL[i].index);
				//~ cout<<vectL[i].index<<end;
				readStructsVec[vectL[i].index].takePref=false;
			}
		}
		++i;
	}
	uint j(0);
	while (j < vectR.size()){
		if (seqsToRemove.unordered_set::count(vectR[j].sequence)){
			if (not vectR[j].canonical ){  // prefix
				//~ readsToRemove.insert(vectR[j].index);
				//~ cout<<vectR[i].index<<endl;
				readStructsVec[vectR[j].index].takePref=false;
			}
		}
		++j;
	}
}


/* from sequences spotted, get reads' index for reads we do not want the suffix to be present in the next pass */
void appendListReadsToRemoveSuff(const vector<edge>& vectL, const vector<edge>& vectR, const unordered_set<string>& seqsToRemove, vector<readStruct>& readStructsVec){
	vector<edge> vectResult;
	uint i(0);
	while (i < vectL.size()){
		if (seqsToRemove.unordered_set::count(vectL[i].sequence)){
			if (vectL[i].canonical == false){  // suffix
				//~ readsToRemove.insert(vectL[i].index);
				readStructsVec[vectL[i].index].takeSuff=false;
			}
		}
		++i;
	}
	uint j(0);
	while (j < vectR.size()){
		if (seqsToRemove.unordered_set::count(vectR[j].sequence)){
			if (vectR[j].canonical == true){  // suffix
				//~ readsToRemove.insert(vectR[j].index);
				readStructsVec[vectR[j].index].takeSuff=false;
			}
		}
		++j;
	}
}



bool compareEdgeByString(const edge& seqL, const edge& seqR){
    return seqL.sequence < seqR.sequence;
}


/*  compaction of two unitigs if they have a k-overlap */
string compaction( readStruct& seq1,  readStruct& seq2, uint k){
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
			//~ if (readsToRemovePref.unordered_set::count(seq1.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq1.takePref){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			//~ if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				//~ readsToRemovePref.insert(seq1.index);
			//~ }
			if (not seq2.takeSuff){
				seq1.takePref=false;
			}else{
				seq1.takePref=true;
			}
			return compactedToReturn;
		} else { // register read's new suffix (the sequence 2's suffix)
			//~ if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq2.takeSuff){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			return compacted;
		}
	} else if (end2.sequence == beg1.sequence) { //  overlap RR
		string compacted(seq2.sequence + seq1.sequence.substr(k));
		bool b(isCanonical(compacted));
		if (not b){
			string compactedToReturn(getCanonical(compacted));
			//~ if (readsToRemoveSuff.unordered_set::count(seq1.index)){
				//~ readsToRemovePref.insert(seq1.index);
			//~ }
			if (not seq1.takeSuff){
				seq1.takePref=false;
			}else{
				seq1.takePref=true;
			}
			//~ if (readsToRemovePref.unordered_set::count(seq2.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq2.takePref){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			return compactedToReturn;
		} else { // register read's new prefix
			return compacted;
		}
	} else if (rEnd1.sequence == end2.sequence) { //  overlap FR
		string compacted(seq1.sequence + rSeq2.substr(k));
		bool b(isCanonical(compacted));
		if (not b){ // new suffix is the read 1's prefix, new prefix is read 2's prefix
			string compactedToReturn(getCanonical(compacted));
			//~ if (readsToRemovePref.unordered_set::count(seq1.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq1.takePref){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			//~ if (readsToRemovePref.unordered_set::count(seq2.index)){
				//~ readsToRemovePref.insert(seq1.index);
			//~ }
			if (not seq2.takePref){
				seq1.takePref=false;
			}else{
				seq1.takePref=true;
			}
			return compactedToReturn;
		} else { // register read's new suffix which is read 2's prefix
			//~ if (readsToRemovePref.unordered_set::count(seq2.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq2.takePref){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			return compacted;
		}
	} else if (beg1.sequence == rBeg2.sequence){ //  overlap RF
		string compacted(rSeq2 + seq1.sequence.substr(k));
		bool b(isCanonical(compacted));
		if (not b){ // new read's suffix is read 2's suffix, new read's prefix is read 1's suffix
			string compactedToReturn(getCanonical(compacted));
			//~ if (readsToRemoveSuff.unordered_set::count(seq1.index)){
				//~ readsToRemovePref.insert(seq1.index);
			//~ }
			if (not seq1.takeSuff){
				seq1.takePref=false;
			}else{
				seq1.takePref=true;
			}
			//~ if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				//~ readsToRemoveSuff.insert(seq1.index);
			//~ }
			if (not seq2.takeSuff){
				seq1.takeSuff=false;
			}else{
				seq1.takeSuff=true;
			}
			return compactedToReturn;
		} else { // register read's new prefix which is read 2's suffix
			//~ if (readsToRemoveSuff.unordered_set::count(seq2.index)){
				//~ readsToRemovePref.insert(seq1.index);
			//~ }
			if (not seq2.takeSuff){
				seq1.takePref=false;
			}else{
				seq1.takePref=true;
			}
			return compacted;
		}
	} else {
		cout<<"fail..."<<endl;
		return "";
	}
}


/*  Recursive func to compact unitigs (if they are sequences or indexes)). If two unitigs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate. */
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k){
	if (not vec[indexreadStruct1].sequence.empty()){
		if (not vec[indexreadStruct2].sequence.empty()){
			string c = compaction(vec[indexreadStruct1], vec[indexreadStruct2], k);
			if (not c.empty()){
				vec[indexreadStruct1].sequence = c;
				vec[indexreadStruct2].index = vec[indexreadStruct1].index;
				vec[indexreadStruct2].sequence = "";
			}
		} else {
			compactInVector(vec, indexreadStruct1, vec[indexreadStruct2].index, k); //  each time a sequence is empty, the index leads to the sequence it's been compacted in-> recursive call until we find the sequence
		}
	} else {
		if (not vec[indexreadStruct2].sequence.empty()){
			compactInVector(vec, vec[indexreadStruct1].index, indexreadStruct2, k);
		} else {
			compactInVector(vec, vec[indexreadStruct1].index, vec[indexreadStruct2].index, k);
		}
	}
}


/*  checks from the suffixes and prefixes of pairs of unitigs if they can be compacted. Appends overlaps that should be removed in the next pass. */
void parseVector(vector<vector<edge>> & left, vector<vector<edge>>& right, vector<readStruct>& readStructsVec, uint k, unordered_set<string>& seqsToRemoveInSuff, unordered_set<string>& seqsToRemoveInPref){
	uint compac(0);
	//SORT
	//~ sort(left.begin(), left.end(), compareEdge());
	//~ sort(right.begin(), right.end(), compareEdge());
	vector<edge> leftSingles;
	vector<edge> rightSingles;
	//~ cout<<"go"<<endl;
	for(uint i(0);i<nBucketsOverlap;++i){
		//~ cout<<left[i].size()<<" "<<right[i].size()<<endl;
		sort(left[i].begin(), left[i].end(), compareEdge());
		sort(right[i].begin(), right[i].end(), compareEdge());
		//remove duplicate
		if (left[i].size()>1){
			leftSingles = removeNotSinglesInLeft(left[i], seqsToRemoveInPref, seqsToRemoveInSuff);
		} else {
			leftSingles = left[i];
		}
		if (right[i].size()>1){
			rightSingles = removeNotSinglesInRight(right[i], seqsToRemoveInPref, seqsToRemoveInSuff);
		} else {
			rightSingles = right[i];
		}
		//~ cout<<"update"<<endl;
		//update set
		appendListReadsToRemovePref(leftSingles, rightSingles, seqsToRemoveInPref, readStructsVec);
		appendListReadsToRemoveSuff(leftSingles, rightSingles, seqsToRemoveInSuff, readStructsVec);
		uint indexL(0),indexR(0);
		//look up for compaction
		//~ cout<<"compaction"<<endl;
		while (indexL < leftSingles.size() and indexR < rightSingles.size()){
			if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
				if (leftSingles[indexL].index != rightSingles[indexR].index){

					compactInVector(readStructsVec, leftSingles[indexL].index, rightSingles[indexR].index, k);
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
	}
}


void parseVectorTip(vector<vector<edge>> & left, vector<vector<edge>>& right, vector<readStruct>& readStructsVec, uint k, uint minSize){
	uint compac(0);
	//SORT
	vector<edge> leftSingles;
	vector<edge> rightSingles;
	for(uint i(0);i<nBucketsOverlap;++i){
		sort(left[i].begin(), left[i].end(), compareEdge());
		sort(right[i].begin(), right[i].end(), compareEdge());
		//remove duplicate
		left[i].erase( unique( left[i].begin(), left[i].end(), compareEdge() ), left[i].end() );
		right[i].erase( unique( right[i].begin(), right[i].end(), compareEdge() ), right[i].end() );
		uint indexL(0),indexR(0);
		//look up for tips
		while (indexL < leftSingles.size() and indexR < rightSingles.size()){
			if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
				++indexL;
				++indexR;
			} else {
				if (leftSingles[indexL].sequence < rightSingles[indexR].sequence){
					++indexL;
					if(leftSingles[indexL].sequence.size()<minSize){
						readStructsVec[leftSingles[indexL].index].sequence="";
					}
				} else {
					++indexR;
					if(rightSingles[indexR].sequence.size()<minSize){
						readStructsVec[rightSingles[indexR].index].sequence="";
					}
				}
			}
		}
	}
}


/* fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs */
void fillPrefVector(vector<vector <edge>>& vecLeft, vector <vector <edge>>& vecRight, const readStruct& seq, uint k, string& rev, string& canonPrefix){
	if (seq.takePref){
		edge prefix = nPrefix(k, seq.index, seq.sequence, true);
		canonPrefix = getStrictCanonical2(prefix.sequence,rev);
		if (not canonPrefix.empty()){ // if is empty, it means the prefix is the rc of itself, we dont we want add it to the vector
			if (prefix.sequence == canonPrefix){
				vecLeft[strHash(canonPrefix)%nBucketsOverlap].push_back({prefix.index, canonPrefix, true});
				//~ cout<<strHash(canonPrefix)%nBucketsOverlap<<endl;
			} else {
				vecRight[strHash(canonPrefix)%nBucketsOverlap].push_back({prefix.index, canonPrefix, false});
			}
		}
	}
}


/* fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs */
void fillSuffVector(vector<vector <edge>>& vecLeft, vector <vector <edge>>& vecRight, const readStruct& seq, uint k, string& rev, string& canonSuffix){
	if (seq.takeSuff){
		edge suffix = nSuffix(k, seq.index, seq.sequence, true);
		canonSuffix = getStrictCanonical2(suffix.sequence,rev);
		if (not canonSuffix.empty()){ // if is empty, it means the suffix is the rc of itself, we dont we want add it to the vector
			if (suffix.sequence == canonSuffix){
				vecRight[strHash(canonSuffix)%nBucketsOverlap].push_back({suffix.index, canonSuffix, true});
			} else {
				vecLeft[strHash(canonSuffix)%nBucketsOverlap].push_back({suffix.index, canonSuffix, false});
			}
		}
	}
}


/* remove duplicates in reads*/
void cleanDuplicatesInreadStructs(vector <readStruct>& vec){
	//~ return;
	//~ cout<<vec.size()<<endl;
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


void cleanDuplicatesInreadStructs2(vector <readStruct>& vec, uint thresholdCleaning ){
	uint i(1);
	string previousSeq(vec[0].sequence);
	uint good(1);
	uint indiceLastSeq(0);
	string temp;
	cout<<vec.size()<<endl;
	while(i<vec.size()){
		temp = vec[i].sequence;
		if (temp == previousSeq){
			vec[i].sequence = "";
			good++;
		}else{
			if(good<thresholdCleaning){
				cout<<indiceLastSeq<<endl;
				vec[indiceLastSeq].sequence="";
			}
			previousSeq = temp;
			indiceLastSeq=i;
			good=1;
		}
		++i;
	}
	if(good< thresholdCleaning){
		vec[indiceLastSeq].sequence="";
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


void printReadStructsIndex(vector <readStruct>& vec,const string& outfileName){
	ofstream out(outfileName);
	for (uint i(0); i<vec.size(); ++i){
		if (not vec[i].sequence.empty()){
			out<<">"+to_string(i)<<endl;
			out<<vec[i].sequence<<endl;
			//~ cin.get();
		}
	}
}


void initVectofreadStructs(vector <readStruct>& vec, const string& sequence){
	if (not sequence.empty()){
		vec.push_back({0, sequence,true,true});
	}
}
