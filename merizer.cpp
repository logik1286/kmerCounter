#include "merizer.h"
#include <string>

using namespace std;
template<> vector<vector<char>> Merizer<vector<char>>::getMerTokens(vector<char> sequence){

	std::vector<vector<char>> rtn;

	if (sequence.size() < kmerWidth)
		return rtn;

	size_t totalKmers = sequence.size() - kmerWidth + 1;
	 rtn.resize(totalKmers);

	auto begin = sequence.begin();

	for (auto i = 0; i < totalKmers; i++) {

		rtn[i] = vector<char>(begin + i, begin + i + kmerWidth);

	}

	return rtn;

}

template<> vector<char> Merizer<vector<char>>::tokenToSequence(vector<char> token) {

	return token;

}
