#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct edge {
  int index;
  string sequence;
};


struct read {
  int index;
  string sequence;
};


edge nPrefix(uint n, int index, const string& sequence){
	return {index, sequence.substr(0, n)};
}



edge nSuffix(uint n, int index, const string& sequence){
	return {index, sequence.substr(sequence.size()-n, n)};
}



string revComp(const string& seq){
	string revCompSeq = "";
	int pos = seq.size()-1;
	char nt;
	do{
		nt = seq[pos];
		switch (nt) {
			case 'A':
				revCompSeq += 'T';
				break;
			case 'T':
				revCompSeq += 'A';
				break;
			case 'C':
				revCompSeq += 'G';
				break;
			case 'G':
				revCompSeq += 'C';
				break;
		}
		--pos;
	} while (pos>=0);
	return revCompSeq;
}



vector<edge> removeDuplicates(vector<edge>& vect){
	vector<edge> vectResult;
	for (uint i(0); i< vect.size(); ++i){
		if (i == vect.size()-1 or vect[i].sequence!=vect[i+1].sequence){
			vectResult.push_back(vect[i]);
		}
	}
	return vectResult;
}



bool compareEdgeByString(const edge& seqL, const edge& seqR)
{
    return seqL.sequence <= seqR.sequence;
}



bool compareReadByString(const read& seqL, const read& seqR)
{
    return seqL.sequence <= seqR.sequence;
}



string getCanonical(const string& seq){
	string revCompSeq = revComp(seq);
	return min(seq, revCompSeq); 
}


//  compaction of two reads if they have a k-overlap
string compaction(const read& seq1, const read& seq2, int k){
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
		return "";
	}
}


//  If two reads are compacted the result replaces one of them, the other becomes "" and keeps the index of its mate.
void compactInVector(vector<read>& vec, int indexRead1, int indexRead2, int k){
	if (not vec[indexRead1].sequence.empty()){
		if (not vec[indexRead2].sequence.empty()){
			string c = compaction(vec[indexRead1], vec[indexRead2], k);
			if (not c.empty()){
				vec[indexRead1] = {vec[indexRead1].index, c};
				vec[indexRead2].index = vec[indexRead1].index;
				vec[indexRead2].sequence = "";
			}
		} else {
			compactInVector(vec, indexRead1, vec[indexRead2].index, k); //  each time a sequence is empty, the index leads to the sequence it's been compacted in-> recursive call until we find the sequence
		}
	} else { 
		if (not vec[indexRead2].sequence.empty()){
			compactInVector(vec, vec[indexRead1].index, indexRead2, k);
		} else {
			compactInVector(vec, vec[indexRead1].index, vec[indexRead2].index, k);
		}
	}
}


//  checks from the suffixes and prefixes of pairs of reads of a vector if they can be compacted
void parseVector(vector<edge>& left, vector<edge>& right, vector<read>& readsVec, int k){
	sort(left.begin(), left.end(), compareEdgeByString);
	sort(right.begin(), right.end(), compareEdgeByString);
	vector<edge> leftSingles = removeDuplicates(left);
	vector<edge> rightSingles = removeDuplicates(right);
	uint indexL(0);
	uint indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			//~ cout << "compaction of " << rightSingles[indexR].sequence<< " from read " << leftSingles[indexL].index <<  " and read " << rightSingles[indexR].index<<endl;
			compactInVector(readsVec, leftSingles[indexL].index, rightSingles[indexR].index, k);
			++ indexL;
			++ indexR;
		} else {
			if (leftSingles[indexL].sequence <= rightSingles[indexR].sequence){
				++indexL;
			} else {
				++indexR;
			}
		}
	}
}


// fill vectors of prefixes and suffixes with canonical k-mers coming from prefixes of reads
void fillPrefVector(vector <edge>& vecLeft, vector <edge>& vecRight, const read& seq, int k){
	edge prefix = nPrefix(k, seq.index, seq.sequence);
	string canonPrefix = getCanonical(prefix.sequence);
	if (prefix.sequence == canonPrefix){
		vecLeft.push_back({prefix.index, canonPrefix});
	} else {
		vecRight.push_back({prefix.index, canonPrefix});
	}
}


// fill vectors of prefixes and suffixes with canonical k-mers coming from suffixes of reads
void fillSuffVector(vector <edge>& vecLeft, vector <edge>& vecRight, const read& seq, int k){
	edge suffix = nSuffix(k, seq.index, seq.sequence);
	string canonSuffix = getCanonical(suffix.sequence);
	if (suffix.sequence == canonSuffix){
		vecRight.push_back({suffix.index, canonSuffix});
	} else {
		vecLeft.push_back({suffix.index, canonSuffix});
	}
}



void cleanDuplicatesInReads(vector <read>& vec){
	uint i(0);
	string previousSeq("");
	while(i<vec.size()){
		string temp = vec[i].sequence;
		if (temp == previousSeq){
			vec[i].sequence = "";
		}
		previousSeq = temp;
		++i;
	}
}


void setReadsIndex(vector <read>& vec){
	for (uint i(0); i<vec.size(); ++i){
		if (not vec[i].sequence.empty()){
			vec[i].index = i;
		}
	}
}



void initVectofReads(vector <read>& vec, string sequence){
	if (not sequence.empty()){
		vec.push_back({0, sequence});
	}
}



int main(int argc, char ** argv){
	if (argc < 3){
		cout << "command line: ./kMILL reads.fasta k" << endl;
	} else {
		string fileName = argv[1];
		uint k = stoi(argv[2]);
		ifstream readFile(fileName);
		string sequence;
		vector <read> sequencesVec;
		while (not readFile.eof()){
			getline(readFile, sequence);
			getline(readFile, sequence);
			initVectofReads(sequencesVec, sequence);
		}
		sort(sequencesVec.begin(), sequencesVec.end(), compareReadByString);
		cleanDuplicatesInReads(sequencesVec);
		setReadsIndex(sequencesVec);
		do {
			vector <edge> right;  // vector of canonical suffixes
			vector <edge> left; //  vector of canonical prefixes
			for (uint i(0); i<sequencesVec.size(); ++i){
				if (sequencesVec[i].sequence.size() > k){
					fillPrefVector(left, right, sequencesVec[i], k);
					fillSuffVector(left, right, sequencesVec[i], k);
				}
			}
			parseVector(left, right, sequencesVec, k);
			--k;
		} while (k>2); 
		for (uint i(0); i<sequencesVec.size(); ++i){
			if (not sequencesVec[i].sequence.empty()){
				cout << sequencesVec[i].index <<  " " <<sequencesVec[i].sequence << endl;
			}
		}
	}
	return 0;
}
