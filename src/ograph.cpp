#include "ograph.h"
#include "utils.h"
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
	if (iL != iR){
		//~ cout<<"go compaction !!!!!!!!!!!!!!!!"<<endl;
		uint s1(unitigs[iL].size()),s2(unitigs[iR].size());
		bool b1(isNumber(unitigs[iL][0])),b2(isNumber(unitigs[iR][0]));
		if(b1 and b2){return compaction(stoi(unitigs[iL]),stoi(unitigs[iR]));}
		if(b1){return compaction(stoi(unitigs[iL]),iR);}
		if(b2){return compaction(iL,stoi(unitigs[iR]));}

		kmer beg1(beg2int128(unitigs[iL]));
		kmer end2(end2int128(unitigs[iR]));

		if(beg1==end2){
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer endrc2(beg2int128rc(unitigs[iR]));
		if(beg1==endrc2){
			reverseinplace2(unitigs[iR]);
			unitigs[iR]+=(unitigs[iL].substr(k));
			unitigs[iL]=to_string(iR);
			return;
		}

		kmer beg2(rcb(endrc2));
		kmer end1(end2int128(unitigs[iL]));
		if(end1==beg2){
			unitigs[iL]+=(unitigs[iR].substr(k));
			unitigs[iR]=to_string(iL);
			return;
		}

		kmer begrc2(rcb(end2));
		if(end1==begrc2){
			unitigs[iL]+=(reverseinplace(unitigs[iR]).substr(k));
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

		kL=left[iL];
		kR=right[iR];
		//~ if(getCanonical(unitigs[kL.indice])==getCanonical("AAAATAAGCCAATACGATCTCAACGCTATTGAAGCGGCTTGCCAGCTAAAGCAACAGGCAGCAGAGGCGCAGGTGACAGCCTTAAGTGTGGGCGGTAAAG")){
			//~ cout<<"hey"<<endl;
		//~ }
		//~ if(getCanonical(unitigs[kR.indice])==getCanonical("AAAATAAGCCAATACGATCTCAACGCTATTGAAGCGGCTTGCCAGCTAAAGCAACAGGCAGCAGAGGCGCAGGTGACAGCCTTAAGTGTGGGCGGTAAAG")){
			//~ cout<<"ho"<<endl;
		//~ }
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
			if(go){compaction(kL.indice,kR.indice);}
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
	if(true){
		kmer kmer1(beg2int128(unitigs[indiceUnitigs]));
		kmer kmer2(rcb(kmer1));
		if(kmer1<kmer2){
			left.push_back(kmerIndice{indiceUnitigs,kmer1});
		}else{
			right.push_back(kmerIndice{indiceUnitigs,kmer2});
		}
	}
	if(true){
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

