#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include "utils.h"
#include "compaction.h"
#include "ograph.h"


using namespace std;


edge nPrefix(uint n, uint index, const string& sequence){
	return {index, sequence.substr(0, n)};
}


edge nSuffix(uint n, uint index, const string& sequence){
	return {index, sequence.substr(sequence.size()-n, n)};
}


vector<edge> removeNotSingles(const vector<edge>& vect, uint k){
	vector<edge> vectResult;
	uint i(0);
	while (i < vect.size()){
		if (i == 0){
			if (vect[i].sequence != vect[i+1].sequence){
				vectResult.push_back(vect[i]);
			}
		} else if (i == vect.size()-1){
			if (vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			}
		} else {
			if (vect[i].sequence != vect[i+1].sequence and vect[i].sequence != vect[i-1].sequence){
				vectResult.push_back(vect[i]);
			}
		}
		++i;
	}
	return vectResult;
}


bool compareEdgeByString(const edge& seqL, const edge& seqR){
    return seqL.sequence < seqR.sequence;
}


//  compaction of two readStructs if they have a k-overlap
string compaction(const readStruct& seq1, const readStruct& seq2, uint k){
	edge beg1 = nPrefix(k, seq1.index, seq1.sequence);
	edge beg2 = nPrefix(k, seq2.index, seq2.sequence);
	edge end1 = nSuffix(k, seq1.index, seq1.sequence);
	edge end2 = nSuffix(k, seq2.index, seq2.sequence);
	edge rEnd1 = {end1.index, revComp(end1.sequence)};
	edge rBeg2 = {beg2.index, revComp(beg2.sequence)};
	string rSeq2 = revComp(seq2.sequence);
	if (end1.sequence == beg2.sequence){ //  overlap FF
		return seq1.sequence + seq2.sequence.substr(k);
	} else if (end2.sequence == beg1.sequence) { //  overlap RR
		return seq2.sequence + seq1.sequence.substr(k);
	} else if (rEnd1.sequence == end2.sequence) { //  overlap FR
		return seq1.sequence + rSeq2.substr(k);
	} else if (beg1.sequence == rBeg2.sequence){ //  overlap RF
		return rSeq2 + seq1.sequence.substr(k);
	} else {
		cout<<"fail..."<<endl;
		return "";
	}
}


//  If two readStructs are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate.
void compactInVector(vector<readStruct>& vec, uint indexreadStruct1, uint indexreadStruct2, uint k){
	if (not vec[indexreadStruct1].sequence.empty()){
		if (not vec[indexreadStruct2].sequence.empty()){
			/*db*/
			//~ if (vec[indexreadStruct1].sequence == "AGCACCGGTGATCATGTTTTTAACATAGTCGGCGTGCCCCGGGCAGTCTACGTGTGCGTAGTGACGGGTCGGGGTGTCGTATTCAACGTGAGAAGTGTTG"){
				//~ cout << "*****************************" << vec[indexreadStruct1].index << " " << indexreadStruct2 << " " << k << endl;
			//~ } else if (vec[indexreadStruct2].sequence == "AGCACCGGTGATCATGTTTTTAACATAGTCGGCGTGCCCCGGGCAGTCTACGTGTGCGTAGTGACGGGTCGGGGTGTCGTATTCAACGTGAGAAGTGTTG"){
				//~ cout << "*****************************" << vec[indexreadStruct2].index << " " << indexreadStruct1 << " " << k << endl;
			//~ }
			//~ if ( indexreadStruct1 == 175398 or indexreadStruct2 == 175398 or indexreadStruct1 == 754701 or indexreadStruct2 == 754701){
				//~ cout << indexreadStruct1 << " " << vec[indexreadStruct1].sequence << endl;
				//~ cout << indexreadStruct2 << " " << vec[ indexreadStruct2].sequence << endl;
				//~ }
			/*end*/
			/*db*/
				//~ if(vec[indexreadStruct1].sequence=="AGCACCGGTGATCATGTTTTTAACATAGTCGGCGTGCCCCGGGCAGTCTACGTGTGCGTAGTGACGGGTCGGGGTGTCGTATTCAACGTGAGAAGTGTTG"){
					//~ cout << "seq in read " << indexreadStruct1 << endl;
				//~ } else if(vec[indexreadStruct2].sequence=="AGCACCGGTGATCATGTTTTTAACATAGTCGGCGTGCCCCGGGCAGTCTACGTGTGCGTAGTGACGGGTCGGGGTGTCGTATTCAACGTGAGAAGTGTTG") {
					//~ cout << "seq in read " << indexreadStruct2 << endl;
				//~ }
				if (indexreadStruct1== 4461 or indexreadStruct2== 4461){
					cout << "*****" << endl;
					cout <<  " " <<vec[indexreadStruct1].index << " " << vec[indexreadStruct1].sequence << endl;
					cout <<" " <<vec[indexreadStruct2].index << " "  <<vec[indexreadStruct2].sequence << endl;
					cout << "compaction" << indexreadStruct1  << " " << indexreadStruct2 << endl;
				}
			/*end*/
			string c = compaction(vec[indexreadStruct1], vec[indexreadStruct2], k);
			/* db*/
			string lol = getCanonical(c) ;
			if (lol == "ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG"){
				cout << vec[indexreadStruct1].sequence << endl;
				cout << vec[indexreadStruct2].sequence << endl;
				cout << "c" << endl;
			}

			/*end*/
			if (not c.empty()){

				vec[indexreadStruct1] = {vec[indexreadStruct1].index, c};
				
				vec[indexreadStruct2].index = vec[indexreadStruct1].index;
				vec[indexreadStruct2].sequence = "";
			}
			/*db*/
			if ( indexreadStruct2==4461 ){
				cout <<  "change d'index " <<vec[indexreadStruct1].index << " " << vec[indexreadStruct1].sequence << endl;
				cout <<" " <<vec[indexreadStruct2].index << " "  <<vec[indexreadStruct2].sequence << endl;
			}
			/*end*/
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


//  checks from the suffixes and prefixes of pairs of readStructs of a vector if they can be compacted
void parseVector(vector<edge>& left, vector<edge>& right, vector<readStruct>& readStructsVec, uint k){
	uint compac(0);
	sort(left.begin(), left.end(), compareEdge());
	sort(right.begin(), right.end(), compareEdge());
	vector<edge> leftSingles;
	vector<edge> rightSingles;
	if (left.size()>1){
		leftSingles = removeNotSingles(left,k);
	} else {
		leftSingles = left;
	}
	if (right.size()>1){
		rightSingles = removeNotSingles(right,k);
	} else {
		rightSingles = right;
	}
	uint indexL(0),indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
			/* debug */
			string debug(getCanonical(leftSingles[indexL].sequence));
			//~ cout << "lol " << debug << endl;
			//~ if ( debug ==  "AAATGGATAACTGGATAGTGAAATAATGCGGACACAGTGGCCCTCTCCGGCAAAACTTAATCTGTTTTTATACATTACCGGTCAGCGTGCGGATGGTTA" or rightSingles[indexR].index ==1437){
//~ //>seq_754701
				//~ cout << "compaction: " << k << " " << leftSingles[indexL].index << " " << leftSingles[indexL].sequence  << " " << rightSingles[indexR].index << " " << rightSingles[indexR].sequence << endl;
			//~ }
			/* end */
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			if (leftSingles[indexL].index != rightSingles[indexR].index){

				/*db*/
				/*end*/
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


// fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of readStructs
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k){
	edge prefix = nPrefix(k, seq.index, seq.sequence);
	string canonPrefix = getCanonical(prefix.sequence);
	if (canonPrefix ==  "AAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTG"){
		cout << "pref" << endl;
	}
	//~ cout<<"pre"<<endl;
	//~ cout<<prefix.sequence<<" "<<canonPrefix<<endl;
	if (prefix.sequence == canonPrefix){
		vecLeft.push_back({prefix.index, canonPrefix});
	} else {
		vecRight.push_back({prefix.index, canonPrefix});
	}
}


// fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of readStructs
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const readStruct& seq, uint k){
	edge suffix = nSuffix(k, seq.index, seq.sequence);
	string canonSuffix = getCanonical(suffix.sequence);
	//~ cout<<"suf"<<endl;
	if (canonSuffix == "AAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTG"){
		cout << "suff" << endl;
	}
	//~ cout<<suffix.sequence<<" "<<canonSuffix<<endl;
	if (suffix.sequence == canonSuffix){
		vecRight.push_back({suffix.index, canonSuffix});
	} else {
		vecLeft.push_back({suffix.index, canonSuffix});
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
			//~ cout << vec[i].index << endl;
		//~ } else {
			//~ cout << "empty" << endl;
		}
	}
}


void initVectofreadStructs(vector <readStruct>& vec, const string& sequence){
	if (not sequence.empty()){
		vec.push_back({0, sequence});
	}
}
