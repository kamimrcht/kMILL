#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>



using namespace std;


char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}


bool isCanonical(const string& seq){
	if (seq.size() > 1){
		char first(seq[0]);
		char last(revCompChar(seq[seq.size()-1]));
		if (first < last){
			return true;
		} else {
			if (first == last){
				string seq2(seq.substr(1,seq.size()-2));
				return isCanonical(seq2);
			} else {
				return false;
			}
		}
	} else {
		if (seq.size() == 1){
			switch(seq[0]){
				case 'A': return true;
				case 'C': return true;
			}
			return false;
		} else {
			return true;
		}
	}
}


// optimized
string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}


string reverseComplements(const string& seq){
	string revCompSeq = "";
	int pos = seq.size()-1;
	do{
		switch (seq[pos]) {
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


string getCanonical(const string& seq){
	return min(seq,revComp(seq));
}


/* read generation */

char randNuc(){
	switch (rand()%4){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	}
	return 'A';
}


string mutate(string& read, uint n){
	for(uint i(0); i<n; ++i){
		int position(rand()%read.size());
		read[position]=randNuc();
	}
	return read;
}


void randGenome(uint size){
	ofstream out("simulGenome",ios::trunc);
	string res;
	for(uint i(0);i<size;++i){
		res+=randNuc();
	}
	out << "genome|size" << size << endl;
	out << res << endl;
}

/* random reads */
void createinputlm(uint lr,uint k){
	ofstream out("simul",ios::trunc);
	uint r;
	string c;
	string kmer(k,'a');
	for(uint b(0);b<k;b++){
		r=rand()%4;
		switch(r){
			case 1:{
				kmer[b]='A';
				break;
			}
			case 2:{
				kmer[b]='C';
				break;
			}
			case 3:{
				kmer[b]='G';
				break;
			}
			case 0:{
				kmer[b]='T';
				break;
			}
		}
	}
	for(int64_t b(0);b<lr;b++){
		kmer=kmer.substr(1,k-1);
		r=rand()%4;
		switch(r){
			case 1:{
				c='A';
				break;
			}
			case 2:{
				c='C';
				break;
			}
			case 3:{
				c='G';
				break;
			}
			case 0:{
				c='T';
				break;
			}
		}
		kmer+=c;
		if(rand()%100==0){
			if(kmer<revComp(kmer)){
				out<<">lol"<<endl;
				out<<kmer<<endl;
			}else{
				out<<">lol"<<endl;
				out<<revComp(kmer)<<endl;
			}
		}
	}
}


void perfectsReadsFromRef(const string& refName,uint length,uint nbRead){
	string ref;
	ofstream out("simul",ios::trunc);
	ifstream in(refName);
	getline(in,ref);
	getline(in,ref);
	for(uint i(0);i<nbRead;++i){
		out<<">"<<endl;
		out<<(ref.substr(rand()%(ref.size()-length),length))<<endl;
	}
}


/* debug */
bool isSubSequenceInSequence(const string& subseq, const string& seq){
	bool present = false;
	if (not seq.empty()){
		for (uint i(0); i < seq.size(); ++i){
			uint w(0);
			do{
				string subseqToCheck(getCanonical(seq.substr(w,100)));
				++w;
				string canon(getCanonical(subseq));
				if (canon == subseqToCheck){
					present = true;
				}
			} while(w + 100 < seq.size());
		}
	}
	return present;
}
/*end debug*/
