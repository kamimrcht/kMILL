#include "ograph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <thread>


/*
 * constructs the compressed dBG from a list of arbitrary sequences
 */


using namespace std;


string reverseinplace(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
	return str;
}


void reverseinplace2(string& str){
	uint i(str.size()-1),j(0);
	for(; j<str.size()/2; --i, ++j){
		str[i] ^= 4;
		str[j] ^= 4;
		if ((str[i]&3) != 3){str[i]^= 17;}
		if ((str[j]&3) != 3){str[j]^= 17;}
		swap(str[i],str[j]);
	}
	if(str.size()%2==1){
		str[j] ^= 4;
		if ((str[j]&3) != 3){str[j]^= 17;}
	}
}


uint chartoint(char c){
	char d = (c >> 1) & 3;
	if (d > 1)
		d ^= 1;
	return d;
}
//
//
// uint64_t stringtoint(const string& str){
// 	uint64_t res(0);
// 	for(uint i(0);i<str.size();i++){
// 		res<<=2;
// 		res+=chartoint(str[i]);
// 	}
// 	return res;
// }
//
//
// uint64_t stringtointc(const string& str){
// 	uint64_t res(0);
// 	for(int32_t i(str.size()-1);i>=0;i--){
// 		res<<=2;
// 		res+=3-chartoint(str[i]);
// 	}
// 	return res;
// }
//
//
// kmer stringtoint128(const string& str){
// 	kmer res(0);
// 	for(uint i(0);i<str.size();i++){
// 		res<<=2;
// 		res+=chartoint(str[i]);
// 	}
// 	return res;
// }
//

// kmer stringtointc128(const string& str){
// 	kmer res(0);
// 	for(int32_t i(str.size()-1);i>=0;i--){
// 		res<<=2;
// 		res+=3-chartoint(str[i]);
// 	}
// 	return res;
// }
//
//
// uint64_t string2intmin(const string& str){
// 	return min(stringtoint(str),stringtointc(str));
// }
//
//
// bool accordtomin(int min, int left_or_right_min){
// 	if(min == -1){
// 		return true;
// 	}
//
// 	if(left_or_right_min==min)
// 		return true;
//
// 	return false;
//
// }
//
//
// string compaction2(const string& seq1,const string& seq2, int k){
// 	size_t s1(seq1.size()),s2(seq2.size());
// 	if(s1==0){return seq2;}
// 	if(s2==0){return seq1;}
//
// 	string rc2(reversecompletment(seq2));
// 	string rc1(reversecompletment(seq1));
//
// 	if(seq1.substr(0,k)==seq2.substr(s2-k,k)){
// 		return seq2+seq1.substr(k);
// 	}else{
// 		if(rc2.substr(s2-k,k)==seq1.substr(0,k)){
// 			return rc2+seq1.substr(k);
// 		}
// 	}
//
// 	if(seq2.substr(0,k)==seq1.substr(s1-k,k)){
// 		return seq1+seq2.substr(k);
// 	}else{
// 		if(rc1.substr(s1-k,k)==seq2.substr(0,k)){
// 			return rc1+seq2.substr(k);
// 		}
// 	}
//
// 	if(rc1.substr(0,k)==seq2.substr(s2-k,k)){
// 			return seq2+rc1.substr(k);
// 	}else{
// 		if(rc2.substr(s2-k,k)==rc1.substr(0,k)){
// 			return rc2+rc1.substr(k);
// 		}
// 	}
//
// 	if(rc2.substr(0,k)==seq1.substr(s1-k,k)){
// 		return seq1+rc2.substr(k);
// 	}else{
// 		if(rc1.substr(s1-k,k)==rc2.substr(0,k)){
// 			return rc1+rc2.substr(k);
// 		}
// 	}
// 	return seq1;
// }
//
//
// string compaction(const string& seq1,const string& seq2, int k){
// 	int s1(seq1.size()),s2(seq2.size());
// 	if(s1==0 or s2==0){
// 		return seq1;
// 	}
// 	string beg1(seq1.substr(0,k));
// 	if(beg1==seq2.substr(s2-k,k)){
// 		return seq2+seq1.substr(k);
// 	}
//
// 	string end1(seq1.substr(s1-k,k));
// 	if(seq2.substr(0,k)==end1){
// 		return seq1+seq2.substr(k);
// 	}
//
// 	string rc2(reversecompletment(seq2));
// 	if(rc2.substr(s2-k,k)==beg1){
// 		return rc2+seq1.substr(k);
// 	}
//
// 	if(rc2.substr(0,k)==end1){
// 		return seq1+rc2.substr(k);
// 	}
// 	return seq1;
// }
//
//
// string compactionBeg(const string& seq1,const string& seq2, int k){
// 	int s2(seq2.size());
// 	string rc2(reversecompletment(seq2));
// 	//~ string rc1(reversecompletment(seq1));
// 	string beg(seq1.substr(0,k));
// 	string begRc(reversecompletment(beg));
//
// 	if(beg==seq2.substr(s2-k,k)){
// 		return seq2+seq1.substr(k);
// 	}else{
// 		if(beg==rc2.substr(s2-k,k)){
// 			return rc2+seq1.substr(k);
// 		}
// 	}
//
// 	if(begRc==seq2.substr(0,k)){
// 		return rc2.substr(0,s2-k)+seq1;
// 	}else{
// 		if(begRc==rc2.substr(0,k)){
// 			return seq2.substr(0,s2-k)+seq1;
// 		}
// 	}
// 	return seq1;
// }
//
//
// string compactionEnd(const string& seq1,const string& seq2, int k){
// 	int s1(seq1.size()),s2(seq2.size());
// 	string rc2(reversecompletment(seq2));
// 	string end(seq1.substr(s1-k,k));
// 	string endRc(reversecompletment(end));
//
// 	if(end==seq2.substr(0,k)){
// 		return seq1+seq2.substr(k);
// 	}else{
// 		if(end==rc2.substr(0,k)){
// 			return seq1+rc2.substr(k);
// 		}
// 	}
//
// 	if(endRc==seq2.substr(s2-k,k)){
// 		return seq1+rc2.substr(k);
// 	}else{
// 		if(endRc==rc2.substr(s2-k,k)){
// 			return seq1+seq2.substr(k);
// 		}
// 	}
// 	return seq1;
// }
//
//
bool isNumber(char c){return (c<64);}


kmer graph3::beg2int128(const string& str){
	kmer resBeg(0);
	for(uint i(0);i<k;++i){
		resBeg<<=2;
		resBeg+=chartoint(str[i]);
	}
	return resBeg;
}


kmer graph3::beg2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[i]);
	}
	return res;
}


kmer graph3::end2int128rc(const string& str){
	kmer res(0);
	for(int i(k-1);i>=0;i--){
		res<<=2;
		res+=3-chartoint(str[str.size()-k+i]);
	}
	return res;
}


kmer graph3::end2int128(const string& str){
	kmer resEnd(0);
	for(uint i(0);i<k;++i){
		resEnd<<=2;
		resEnd+=chartoint(str[str.size()-k+i]);
	}
	return resEnd;
}


kmer graph3::rcb(kmer min){
	kmer resrcb(0);
	kmer offsetrcb(1);
	for(uint i(0); i<k;++i){
		resrcb+=(3-(min%4))<<(2*(k-1-i));
		min>>=2;
	}
	return resrcb;
}


void graph3::compaction(uint iL,  uint iR){
	// if(iL==lol and found ){
	// 	cout<<"compactionL !"<<endl;
	// }
	// if(iR==lol  and found){
	// 	cout<<"compactionR !"<<endl;
	// }
	/* db*/
	if (unitigs[iR]== "CAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" or unitigs[iL]=="CAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" ){
		cout << "elela" << endl;
	}
	//~ if (unitigs[iR]== "" or unitigs[iL]==""){
		//~ cout << "elela" << endl;
	//~ }
	/*end*/
	if (iL != iR){
		uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
		bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
		if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
		if(b1){return compaction(stoi(unitigs[iL]),iR);}
		if(b2){return compaction(iL,stoi(unitigs[iR]));}

		kmer beg1(beg2int128(unitigs[iL]));
		kmer end2(end2int128(unitigs[iR]));

		if(beg1==end2){
			/* db*/
			if (unitigs[iR]== "ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG"  or unitigs[iL]=="ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" ){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			if (unitigs[iR]== "CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT" or unitigs[iL]=="CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT"){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			/*end*/
			unitigs[iR]+=(unitigs[iL].substr(k));
			// if((iL==lol or iR==lol) and found){
			// 	cout<<unitigs[iR]<<endl;
			// 	lol=iR;
			// }
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer endrc2(beg2int128rc(unitigs[iR]));
		if(beg1==endrc2){
			/* db*/
			if (unitigs[iR]== "ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG"  or unitigs[iL]=="ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" ){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			if (unitigs[iR]== "CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT" or unitigs[iL]=="CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT"){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			/*end*/
			reverseinplace2(unitigs[iR]);
			unitigs[iR]+=(unitigs[iL].substr(k));
			// 	if((iL==lol or iR==lol) and found){
			// 	cout<<unitigs[iR]<<endl;
			// 	lol=iR;
			// }
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer beg2(rcb(endrc2));
		kmer end1(end2int128(unitigs[iL]));
		if(end1==beg2){
			/* db*/
			if (unitigs[iR]== "ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG"  or unitigs[iL]=="ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" ){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			if (unitigs[iR]== "CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT" or unitigs[iL]=="CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT"){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			/*end*/
			unitigs[iL]+=(unitigs[iR].substr(k));
			// 	if((iL==lol or iR==lol) and found){
			// 	cout<<unitigs[iL]<<endl;
			// 	lol=iL;
			// }
			unitigs[iR]=to_string(iL);
			return;
		}

		kmer begrc2(rcb(end2));
		if(end1==begrc2){
			/* db*/
			if (unitigs[iR]== "ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG"  or unitigs[iL]=="ACAGAGTTGGATCCCGGTCGTTTCTGGATTTTTGTTAAGCCGGGTTATTCGTAAAGATCGCGATGAGCCGTTTGTTTATAACAGCGGCGGTTCCTCTTTTG" ){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			if (unitigs[iR]== "CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT" or unitigs[iL]=="CAAAAGAGGAACCGCCGCTGTTATAAACAAACGGCTCATCGCGATCTTTACGAATAACCCGGCTTAACAAAAATCCAGAAACGACCGGGATCCAACTCTGT"){
				cout << "compaction malfoy " << unitigs[iL] << " " << unitigs[iR] << endl;
				if (unitigs[iL] == unitigs[iR]){
					cout << "BORDEL" << endl;
				}
			}
			/*end*/
			unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
			// if((iL==lol or iR==lol) and found){
			// 	cout<<unitigs[iL]<<endl;
			// 	lol=iL;
			// }
			unitigs[iR]=to_string(iL);
			return;
		}
	}
	// cout<<"wut"<<endl;
}


void graph3::debruijn(){
	sort(left.begin(),left.end(),comparator());
	sort(right.begin(),right.end(),comparator());
	uint iL(0),iR(0);
	kmerIndice kL,kR;
	while(iL!=left.size() and iR!=right.size()){
		//~ if (unitigs[iL] == "AAAATGGATAACTGGATAGTGAAATAATGCGGACACAGTGGCCCTCTCCGGCAAAACTTAATCTGTTTTTATACATTACCGGTCAGCGTGCGGATGGTTA" or unitigs[iR] == "AAAATGGATAACTGGATAGTGAAATAATGCGGACACAGTGGCCCTCTCCGGCAAAACTTAATCTGTTTTTATACATTACCGGTCAGCGTGCGGATGGTTA" or  unitigs[iL] == "TAACCATCCGCACGCTGACCGGTAATGTATAAAAACAGATTAAGTTTTGCCGGAGAGGGCCACTGTGTCCGCATTATTTCACTATCCAGTTATCCATTTT"  or  unitigs[iR] ==  "TAACCATCCGCACGCTGACCGGTAATGTATAAAAACAGATTAAGTTTTGCCGGAGAGGGCCACTGTGTCCGCATTATTTCACTATCCAGTTATCCATTTT"){
			//~ cout << "yes" << endl;
			//~ cin.get();
			
		//~ }
		kL=left[iL];
		kR=right[iR];
		if(kL.kmmer==kR.kmmer){
			bool go(true);
			++iL;++iR;
				if(left[iL].kmmer==kL.kmmer){
					go=false;
					while(left[++iL].kmmer==kL.kmmer){}
				}
				if(right[iR].kmmer==kL.kmmer){
					go=false;
					while(right[++iR].kmmer==kR.kmmer){}
				}
			if(go and kL.indice != kR.indice){compaction(kL.indice,kR.indice);}
		}else{
			if(kL.kmmer<kR.kmmer){
				while(left[++iL].kmmer==kL.kmmer){}
			}else{
				while(right[++iR].kmmer==kR.kmmer){}
			}
		}
	}
}


bool graph3::output(uint i){return !isNumber(unitigs[i][0]);}


bool graph3::clear(){delete [] unitigs;return true;}


uint graph3::size(){return indiceUnitigs;};


void graph3::addtuple(tuple<string,uint,uint> tuple){
	unitigs[indiceUnitigs]=move(get<0>(tuple));
	if(minimizer==(get<1>(tuple))){
		kmer kmer1(beg2int128(unitigs[indiceUnitigs]));
		kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			left.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	if(minimizer==get<2>(tuple)){
		kmer kmer1(end2int128(unitigs[indiceUnitigs]));
		kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			right.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			left.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	++indiceUnitigs;
}


// void compareUnitigs(const string& fileFa,const string& fileDot){
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq;
// 	getline(streamFa,seq);
// 	while (!streamFa.eof()) {
// 		getline(streamFa,seq,'>');
// 		seq=seq.substr(0,seq.size()-1);
// 		setFa.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		getline(streamFa,seq);
// 		++c;
// 	}
// 	cout<<1<<endl;
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 		seq=seq.substr(0,seq.size()-1);
// 		setDot.insert(seq);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		++d;
// 	}
// 	cout<<2<<endl;
// 	for(auto it(setFa.begin());it!=setFa.end();++it){
// 		if(setDot.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<3<<endl;
// 	for(auto it(setDot.begin());it!=setDot.end();++it){
// 		if(setFa.count(*it)==0){
// 			++a;
// 		}
// 	}
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }
//
//
// void compareKmers(const string& fileFa,const string& fileDot){
// 	uint k(31);
// 	string kmer;
// 	uint a(0),b(0),c(0),d(0);
// 	unordered_set<string> setFa,setDot;
// 	ifstream streamFa(fileFa),streamDot(fileDot);
// 	string seq,inter,nimp;
//
//
//
// 	// cout<<1<<endl;
// 	while (!streamFa.eof()) {
// 		getline(streamFa,nimp);
// 		// cout<<"nimp"<<nimp<<endl;
// 		getline(streamFa,seq);
// 		// cout<<"seq"<<seq<<endl;
// 		point:
// 		char c=streamFa.peek();
// 		if(c=='>'){
// 			point2:
// 			// seq=seq.substr(0,seq.size());
// 			// for(uint j(0);(j)<seq.size();++j){
// 			// 	if(seq[j]!='A' and seq[j]!='C' and seq[j]!='T' and seq[j]!='G'){
// 			// 		cout<<seq<<endl;
// 			// 		cout<<"lol"<<endl;
// 			// 		exit(0);
// 			// 	}
// 			// }
// 			for (uint i = 0; i+k <=seq.size(); ++i) {
// 				kmer=seq.substr(i,k);
// 				// cout<<kmer<<endl;
// 				kmer=getRepresent(kmer);
// 				// if(setDot.count(kmer)==0){
// 				// 	++a;
// 				// }
// 				setFa.insert(kmer);
// 			}
// 		}else{
// 			if(!streamFa.eof()){
// 				// cout<<"inter"<<endl;
// 				// cout<<seq<<endl;
// 				getline(streamFa,inter);
// 				// cout<<inter<<endl;
// 				seq+=inter;
// 				goto point;
// 			}else{
// 				// cout<<"lol2"<<endl;
// 				goto point2;
// 			}
// 		}
// 	}
// 	cout<<2<<endl;
//
// 	while (!streamDot.eof()){
// 		getline(streamDot,seq);
// 		seq=seq.substr(0,k);
// 		// cout<<seq<<endl;
// 		// cin.get();
// 		if(setFa.count(getRepresent(seq))==0){
// 			cout<<seq<<endl;
// 			++a;
// 		}
// 	}
//
// 	// while (!streamDot.eof()){
// 	// 	getline(streamDot,seq);
// 	// 	transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
// 	// 	seq=seq.substr(0,seq.size()-1);
// 	// 	// cout<<seq<<endl;
// 	// 	for (uint i = 0; i+k <=seq.size(); ++i) {
// 	// 		kmer=seq.substr(i,k);
// 	// 		// cout<<kmer<<endl;
// 	// 		kmer=getRepresent(kmer);
// 	// 		// setDot.insert(kmer);
// 	// 		if(setFa.count(kmer)==0){
// 	// 			++b;
// 	// 		}
// 	// 	}
// 	// 	// cout<<seq<<endl;
// 	// 	// cin.get();
// 	// 	// ++d;
// 	// }
// 	// for(auto it(setFa.begin());it!=setFa.end();++it){
// 	// 	if(setDot.count(*it)==0){
// 	// 		++a;
// 	// 	}
// 	// }
// 	cout<<3<<endl;
// 	// for(auto it(setDot.begin());it!=setDot.end();++it){
// 	// 	if(setFa.count(*it)==0){
// 	// 		++b;
// 	// 	}
// 	// }
// 	cout<<a<<" "<<b<<endl;
// 	cout<<c<<" "<<d<<endl;
// }
