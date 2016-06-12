#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "compaction.h"


using namespace std;


struct pair_hash {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};


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


string getStrictCanonical(const string& seq){ // we dont want to compact sequences whose k prefix or suffix is the rev comp of itself (because of repeats)
	string r(revComp(seq));
	if(seq != r){
		return min(seq, r);
	} else {
		return "";
	} 
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
	for(uint i(0); i < n; ++i){
		int position(rand()%read.size());
		read[position] = randNuc();
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

void mutateReadsFromRef(const string& refName,uint length,uint nbRead){
	string ref;
	ofstream out("simul",ios::trunc);
	ifstream in(refName);
	getline(in,ref);
	getline(in,ref);
	for(uint i(0);i<nbRead;++i){
		out<< ">" <<endl;
		string read(ref.substr(rand()%(ref.size()-length),length));
		if (i%10 == 0){
			read = mutate(read, 1);
		}
		out<< read <<endl;
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


string compactionType(const string& seq1, const string& seq2, uint k){
	edge beg1 = nPrefix(k, 0, seq1, true);
	edge beg2 = nPrefix(k, 0, seq2, true);
	edge end1 = nSuffix(k, 0, seq1, true);
	edge end2 = nSuffix(k, 0, seq2, true);
	edge rEnd1 = {0, revComp(end1.sequence), true};
	edge rBeg2 = {0, revComp(beg2.sequence), true};
	//~ string rSeq2 = revComp(seq2.sequence);
	if (end1.sequence == beg2.sequence){ //  overlap FF
		return "";
	} else if (end2.sequence == beg1.sequence) { //  overlap RR
		return "RR";
	} else if (rEnd1.sequence == end2.sequence) { //  overlap FR
		return "pp";
	} else if (beg1.sequence == rBeg2.sequence){ //  overlap RF
		return "ff";
	} else {
		cout<<"fail..."<<endl;
		return "";
	}
}


void sequences2dot(vector<readStruct>& seqV, uint k, unordered_set<uint>& colorNodePref, unordered_set<uint>& colorNodeSuff, unordered_map<uint, uint>& sizesNode){
    ofstream out("out.dot",ofstream::out);
    out<<"digraph ham {"<<endl;
    // title
    out << "labelloc=\"t\"" << endl ;
    out << "label = \"k="  << k << "\"" <<  endl;
    unordered_multimap<string, readStruct> right2seq;
    unordered_multimap<string, readStruct> left2seq;
    
    string begin, end, sequence;
    
    for(uint i(0);i<seqV.size();++i){
		sequence = seqV[i].sequence;
		if (not sequence.empty()){
			begin = sequence.substr(0,k);
			end = sequence.substr(sequence.size() - k,k);
			if(begin == getCanonical(begin)){
				left2seq.insert({begin, seqV[i]});
			}else{
				right2seq.insert({getCanonical(begin),seqV[i]});
			}
			if(end == getCanonical(end)){

				right2seq.insert({end,seqV[i]});
			}else{
				left2seq.insert({getCanonical(end),seqV[i]});
			}
		}
    }
    //~ for (auto i = left2seq.begin(); i != left2seq.end(); ++i){
		//~ cout << "L " << i->first << endl;
	//~ }
	//~ for (auto i = right2seq.begin(); i != right2seq.end(); ++i){
		//~ cout  << "R "<< i->first << endl;
	//~ }
    uint cbeg;
    uint cend;
     unordered_set<pair<uint,uint>,pair_hash> nadine;
    for(uint i(0);i<seqV.size();++i){
        sequence = seqV[i].sequence;
        begin = sequence.substr(0,k);
        float width(0.005);
        float height(0.005);
        uint notusebeginp(0);
        uint notusebeginf(0);
        uint notuseendf(0);
        uint notuseendp(0);
       
        
        if (not sequence.empty()){
			width *= seqV[i].sequence.size();
			height *= seqV[i].sequence.size();
			if (k>2){
				if(begin == getCanonical(begin)){ // begin would be put in left -> look for overlaps in right
					auto range = right2seq.equal_range(begin);
					for(auto it(range.first);it!=range.second;++it){ // type RR or RF overlaps, node with overlap in R overlaps node with overlap in L
						if (seqV[i].index != it->second.index){
							string type(compactionType(sequence,it->second.sequence, k));
							if (type == "ff" or type == "RR"){
								++ notusebeginf;
							} else {
								++ notusebeginp;
								}
							if(nadine.count({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)})==0){
								if (type != "RR"){
									out << seqV[i].index << "->" << it->second.index << "[ label=\"" << type << "\" ];" << endl;
									nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								} else {
									//~ out << it->second.index << "->" << seqV[i].index << "[ label=\"" << type << "\" ];" << endl;
									//~ nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								}
							}
						}
					}
				}else{ // begin would be in R, look fot overlaps in L
					auto range = left2seq.equal_range(getCanonical(begin));
					for(auto it(range.first);it!=range.second;++it){ // type FF or RF, 1 -> 2
						if (seqV[i].index != it->second.index){
							string type(compactionType(sequence,it->second.sequence, k));
							if (type == "ff" or type == "RR"){
								++ notusebeginf;
							} else {
								++ notusebeginp;
								}
							if(nadine.count({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)})==0){
								string type(compactionType(sequence,it->second.sequence, k));
								if (type != "RR"){
									out << seqV[i].index << "->" << it->second.index << "[ label=\"" << type << "\" ];" << endl;
									nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								} else {
									//~ out << it->second.index << "->" << seqV[i].index << "[ label=\"" << type << "\" ];" << endl;
									//~ ++ cbeg;
									//~ nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								}
							}
						}
					}
				}
				end = sequence.substr(sequence.size()-k,k);
				if(end == getCanonical(end)){ // end would be in R
					auto range = left2seq.equal_range(end);
					for(auto it(range.first);it!=range.second;++it){ // FF or FR, 1-> 2
						if (seqV[i].index != it->second.index){
							string type(compactionType(sequence,it->second.sequence, k));
							if (type == "ff" or type == "RR"){
								++ notuseendf;
							} else {
								++ notuseendp;
								}
							if(nadine.count({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)})==0){
								string type(compactionType(sequence,it->second.sequence, k));
								if (type != "RR"){
									out << seqV[i].index << "->" << it->second.index << "[ label=\"" << type << "\" ];" << endl;
									nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								} else {
									//~ out << it->second.index << "->" << seqV[i].index << "[ label=\"" << type << "\" ];" << endl;
									//~ ++ cbeg;
									//~ nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								}
							}
						}
					}
				}else{
					auto range = right2seq.equal_range(getCanonical(end));
					for(auto it(range.first);it!=range.second;++it){
						if (seqV[i].index != it->second.index){ // RR or FR, 2->1
							string type(compactionType(sequence,it->second.sequence, k));
							if (type == "ff" or type == "RR"){
								++ notuseendf;
							} else {
								++ notuseendp;
								}
							if(nadine.count({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)})==0){
								string type(compactionType(sequence,it->second.sequence, k));
								if (type != "RR"){
									out << seqV[i].index << "->" << it->second.index << "[ label=\"" << type << "\" ];" << endl;
									nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								} else {
									//~ out << it->second.index << "->" << seqV[i].index << "[ label=\"" << type << "\" ];" << endl;
									//~ ++ cbeg;
									//~ nadine.insert({min(seqV[i].index ,it->second.index),max(seqV[i].index ,it->second.index)});
								}
							}
						}
					}
				}
			}
			bool changedSize(false);
			if (sizesNode.count(seqV[i].index)){
				if (sizesNode[seqV[i].index] != seqV[i].sequence.size()){
					changedSize = true;
					sizesNode[seqV[i].index] = seqV[i].sequence.size();
				}
			} else {
				sizesNode.insert({seqV[i].index, seqV[i].sequence.size()});
			}
			if (notusebeginf > 1 or notusebeginp > 1){
				colorNodePref.insert(seqV[i].index);
			}
			if (notuseendf > 1 or notuseendp > 1){
				colorNodeSuff.insert(seqV[i].index);
			}
			if (not changedSize){
				if (colorNodePref.count(seqV[i].index)){
					out << seqV[i].index << "[fillcolor=\"white;0.5:red\" style=filled width=" << width <<  " height=" << height <<  "]" <<endl; // no more compact wit prefix
				}else if (colorNodeSuff.count(seqV[i].index)) {
					out << seqV[i].index << "[fillcolor=\"red;0.5:white\" style=filled width=" << width <<  " height=" << height <<  "]" <<endl; // no more compact with suffix
				}else{
					out << seqV[i].index << "[width=" << width <<  " height=" << height <<  "]" <<endl;
				}
			} else {
				//~ cout << "changed" << endl;
				if (colorNodePref.count(seqV[i].index)){
					out << seqV[i].index << "[shape=rect fillcolor=\"white;0.5:red\" style=filled width=" << width <<  " height=" << height <<  "]" <<endl; // no more compact wit prefix
				}else if (colorNodeSuff.count(seqV[i].index)) {
					out << seqV[i].index << "[shape=rect fillcolor=\"red;0.5:white\" style=filled width=" << width <<  " height=" << height <<  "]" <<endl; // no more compact with suffix
				}else{
					out << seqV[i].index << "[shape=rect width=" << width <<  " height=" << height <<  "]" <<endl;
				}
			}
		}
    }
    out<<"}"<<endl;
}




/*end debug*/
