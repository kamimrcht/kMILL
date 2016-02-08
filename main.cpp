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



//~ string compaction(int index1, int index2, const string& seq1,const string& seq2, int k){
string compaction(const read& seq1, const read& seq2, int k){
	edge beg1 = nPrefix(k, seq1.index, seq1.sequence);
	edge beg2 = nPrefix(k, seq2.index, seq2.sequence);
	edge end1 = nSuffix(k, seq1.index, seq1.sequence);
	edge end2 = nSuffix(k, seq2.index, seq2.sequence);
	edge rEnd1 = {end1.index, revComp(end1.sequence)};
	edge rBeg2 = {beg2.index, revComp(beg2.sequence)};
	string rSeq2 = revComp(seq2.sequence);
	if (end1.sequence == beg2.sequence){ // FF
		return seq1.sequence + seq2.sequence.substr(k);
	} else if (end2.sequence == beg1.sequence) { //RR
		//~ cout << "ce cas" << endl;
		return seq2.sequence + seq1.sequence.substr(k);
	} else if (rEnd1.sequence == end2.sequence) {//FR
		return seq1.sequence + rSeq2.substr(k);
	} else if (beg1.sequence == rBeg2.sequence){//RF
		return rSeq2 + seq1.sequence.substr(k);
	} else {
		return "";
	}
}



void compactInVector(vector<read>& vec, int indexRead1, int indexRead2, int k){
	//~ cout << "compact" << endl;
	if (not vec[indexRead1].sequence.empty()){
		//~ cout << "bli" << endl;
		if (not vec[indexRead2].sequence.empty()){
			//~ cout << vec[indexRead1].sequence << " " << vec[indexRead2].sequence << " " << k << endl;
			string c = compaction(vec[indexRead1], vec[indexRead2], k);
			//~ cout << "bli" << endl;
			if (not c.empty()){
				//~ cout << "bli";
				vec[indexRead1] = {vec[indexRead1].index, c};
				vec[indexRead2].index = vec[indexRead1].index;
				vec[indexRead2].sequence = "";
				//~ cout << vec[indexRead1].index  << vec[indexRead1].sequence << endl;
			}
		} else {
			compactInVector(vec, indexRead1, vec[indexRead2].index, k);
		}
	} else { 
		if (not vec[indexRead2].sequence.empty()){
			compactInVector(vec, vec[indexRead1].index, indexRead2, k);
		} else {
			compactInVector(vec, vec[indexRead1].index, vec[indexRead2].index, k);
		}
	}
}



void parseVector(vector<edge>& left, vector<edge>& right, vector<read>& readsVec, int k){
	sort(left.begin(), left.end(), compareEdgeByString);
	sort(right.begin(), right.end(), compareEdgeByString);
	vector<edge> leftSingles = removeDuplicates(left);
	vector<edge> rightSingles = removeDuplicates(right);
	//~ cout << "***** right removed****" << endl;
	//~ for (uint i(0); i< rightSingles.size();++i){
		//~ cout << rightSingles[i].index << rightSingles[i].sequence << endl;
	//~ }
	//~ cout << "***** left removed****" << endl;
	//~ for (uint i(0); i< leftSingles.size();++i){
		//~ cout << leftSingles[i].index << leftSingles[i].sequence << endl;
	//~ }
	uint indexL(0);
	uint indexR(0);
	while (indexL < leftSingles.size() and indexR < rightSingles.size()){
		//~ cout << indexL << "*" << indexR << endl; 
		//~ cout << leftSingles[indexL].sequence << "-" << rightSingles[indexR].sequence << endl;
		//~ cout << leftSingles[indexL].index << "--" << rightSingles[indexR].index << endl;
		//~ cout << readsVec[leftSingles[indexL].index].index << "--" << readsVec[rightSingles[indexR].index].index << endl;
		//~ cout << readsVec[leftSingles[indexL].index].sequence << "--" << readsVec[rightSingles[indexR].index].sequence << endl;
		if (leftSingles[indexL].sequence == rightSingles[indexR].sequence){
			cout << "compaction of " << rightSingles[indexR].sequence<< " from read " << leftSingles[indexL].index <<  " and read " << rightSingles[indexR].index<<endl;
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
		//~ cout << indexL <<" "<< indexR << endl;
	}
}



int main(int argc, char ** argv){
	ifstream readFile("read.fa");
	string sequence;
	vector <read> sequencesVec;
	//~ int indexCount(0);
	while (not readFile.eof()){
		getline(readFile, sequence);
		getline(readFile, sequence);
		if (not sequence.empty()){
			sequencesVec.push_back({0, sequence});
			//~ ++ indexCount;
		}
	}
	sort(sequencesVec.begin(), sequencesVec.end(), compareReadByString);
	uint i(0);
	string previousSeq("");
	while(i<sequencesVec.size()){
		string temp = sequencesVec[i].sequence;
		if (temp == previousSeq){
			sequencesVec[i].sequence = "";
		}
		previousSeq = temp;
		++i;
	}
	//~ cout << "%%%%%%%%%" << endl;
	for (uint i(0); i<sequencesVec.size(); ++i){
		if (not sequencesVec[i].sequence.empty()){
			sequencesVec[i].index = i;
			//~ cout <<i << " "<< sequencesVec[i].index << sequencesVec[i].sequence << endl;
		}
	}
	//~ cout << "%%%%%%%" <<endl;
	//~ cout << nPrefix(4, sequences[0]) << endl;
	//~ cout << nSuffix(5, sequences[0]) << endl;
	//~ cout << "----" <<endl;
	//~ cout << sequences[0] << ' ' << revComp(sequences[0]) << ' ' << getCanonical(sequences[0]) << endl;

	uint k = 5;
	do {
		vector <edge> right;
		vector <edge> left;
		for (uint i(0); i<sequencesVec.size(); ++i){
			if (sequencesVec[i].sequence.size() > k){
				edge prefix = nPrefix(k, sequencesVec[i].index, sequencesVec[i].sequence);
				edge suffix = nSuffix(k, sequencesVec[i].index, sequencesVec[i].sequence);
				string canonPrefix = getCanonical(prefix.sequence);
				string canonSuffix = getCanonical(suffix.sequence);
				if (prefix.sequence == canonPrefix){
					left.push_back({prefix.index, canonPrefix});
				} else {
					right.push_back({prefix.index, canonPrefix});
				}
				if (suffix.sequence == canonSuffix){
					right.push_back({suffix.index, canonSuffix});
				} else {
					left.push_back({suffix.index, canonSuffix});
				}
			}
		}
		//~ cout << "***** left****" << endl;
		//~ for (uint i(0); i< left.size();++i){
			//~ cout << left[i].index << left[i].sequence << endl;
		//~ }
		//~ cout << "***** right****" << endl;
		//~ for (uint i(0); i< right.size();++i){
			//~ cout << right[i].index << right[i].sequence << endl;
		//~ }
		//~ cout << "****" << endl;
		parseVector(left, right, sequencesVec, k);
		--k;
	} while (k>2); 
	for (uint i(0); i<sequencesVec.size(); ++i){
		if (not sequencesVec[i].sequence.empty()){
			cout << sequencesVec[i].index <<  " " <<sequencesVec[i].sequence << endl;
		}
	}
	
	//~ vector<edge> vL = {{0,"AAA"},{1,"AAA"},{2,"AAC"},{3,"AAT"},{4,"CCC"},{5,"CCC"},{6,"CCC"},{7,"TTG"},{8,"TTT"},{9,"TTT"}};
	//~ vector<edge> vR = {{0,"AAA"},{1,"TTT"}};
	//~ vector<edge> vv = removeDuplicates(vL);
	//~ for (int i(0); i<vv.size(); ++i){
		//~ cout << vv[i].sequence <<endl;
	//~ }
	//~ cout << "----" <<endl;
	//~ parseVector(vL,vR);
	//~ cout << "----" <<endl;

	//~ read f1 = {1, "AATATC"};
	//~ read f2 = {2, "ATCTTT"};
	//~ read r1 = {1, revComp(f1.sequence)}; // GATATT
	//~ read r3 = {3, "ATTGAC"};
	//~ read r4 = {4, "ATCCCCC"};
	//~ read f3 = {5, "ATTACT"};

	//~ read hop = {0, "AATATC"};
	//~ read lol = {0, "CAATAT"};

	//~ string c = compaction(f1, f2, 3);
	//~ string c = compaction(r1, r3, 3);
	//~ string c = compaction(f1, r4, 3);
	//~ string c = compaction(r1, f3, 3);
	//~ string c = compaction(hop, lol, 5);
	//~ string c = compaction({1, "AAAAC"}, {2,"CCCCCC"}, 3);
	//~ cout << c << endl;
	return 0;
}
