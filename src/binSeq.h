#ifndef BINSEQ
#define BINSEQ



#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>



using namespace std;



class binSeq{
public:
	vector<bool> bits;
	binSeq substr(uint begin);
	binSeq substr(uint begin, uint length);
	uint size();
	binSeq(){};
	binSeq(const string& str){
		for(uint i(0);i<str.size();++i){
			switch (str[i])  {
				case 'A':
					bits.push_back(false);
					bits.push_back(false);
				break;
				case 'C':
					bits.push_back(false);
					bits.push_back(true);
				break;
				case 'G':
					bits.push_back(true);
					bits.push_back(false);
				break;
				case 'T':
					bits.push_back(true);
					bits.push_back(true);
				break;
				default:
				break;
			}
		}
	};
	binSeq(const binSeq& bs){
		bits=bs.bits;
	};
};



#endif /* defined(BINSEQ) */
