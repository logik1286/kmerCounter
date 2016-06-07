#pragma once

#include "kmerError.h"
#include <vector>
#include <bitset>
#include <array>

//Base Merizer template. Specialized from this to implement functions, kind of like an interface.
//The Merizer is intended to separate the process of calculating mers from how mers are managed in hash tables or other algorithms,
//so the two can be interchanged without having a lot of duplicate code. 
template <class T> class Merizer {

public:
	Merizer(size_t _kmerWidth) :kmerWidth(_kmerWidth) {};

	//Convert an encoded sequence into kmer tokens. These tokens can be the kmers themselves or a compression of the kmer string 
	std::vector<T> getMerTokens(std::vector<char> encodedSequence);

	//Convert the token into the kmer string
	std::vector<char> tokenToSequence(T token);

private:

	size_t kmerWidth;

};

//http://stackoverflow.com/questions/21245139/fastest-way-to-compare-bitsets-operator-on-bitsets
template<size_t N> struct BitSetComparer {
	bool operator() (const std::bitset<N>& x, const std::bitset<N>& y) const {
		for (int i = N - 1; i >= 0; i--) {
			if (x[i] ^ y[i]) return y[i];
		}
		return false;
	}
};

//This Merizer attempted to leverage that only 3 bits are needed to express GATCN characters,
//and perhaps reduce time spent performing hashing/sorting calculations. It turns out that
//the template library doesn't offer the level of access to really use this compress efficiently,
//so below is another version of this without using STL.
//
//The general idea is that sequential kmers are translations of previous ones, so why not use
//bitwise shifting to leverage this behavior towards fewer computations.
template <size_t N> class Merizer<std::bitset<N>> {

public:

	static size_t const bitsPerBasePair = 3;

	static size_t neededBitsForKmerWidth(size_t kmerWidth) {
		return kmerWidth * bitsPerBasePair;
	};

	Merizer(size_t _kmerWidth) :kmerWidth(_kmerWidth){
	
		bitsNeeded = neededBitsForKmerWidth(kmerWidth);

		if (bitsNeeded > N)
			throw KmerError(1, "Invalid bitset size for provided kmer width. Need at least " + std::to_string(bitsNeeded) + " bits");

		//create mask
		mask.reset();
		for (auto i = 0; i < bitsNeeded; i++) {
			mask[i] = true;
		}

		firstBasePairMask.reset();
		for (auto i = 0; i < bitsPerBasePair; i++) {
			firstBasePairMask[i] = true;
		}

	};

	std::vector<std::bitset<N>> getMerTokens(std::vector<char> sequence) {

		std::vector<std::bitset<N>> rtn;

		if (sequence.size() < kmerWidth)
			return rtn;

		size_t totalKmers = sequence.size() - kmerWidth + 1;
		 rtn.resize(totalKmers);

		std::bitset<N> tokenizer; tokenizer.reset();

		//Need to prime tokenizer
		for (auto i = 0; i < kmerWidth-1; i++)
			pushBasePairToTokenizer(tokenizer, sequence[i]);

		//Now start collecting tokens;
		for (auto i = 0; i < totalKmers; i++) {
			pushBasePairToTokenizer(tokenizer, sequence[i+kmerWidth-1]);
			rtn[i] = tokenizer;
		}

		return rtn;

	};


	std::vector<char> tokenToSequence(std::bitset<N> token) {

		std::vector<char> rtn; rtn.resize(kmerWidth);

		for (auto i = 0; i < kmerWidth; i++) {

			auto thisElement = kmerWidth - 1 - i;
			auto thisChar = token & firstBasePairMask;
			rtn[thisElement] = thisChar;// getCharacterFromValue(thisChar);
			token >>= bitsPerBasePair;

		}

		return rtn;

	};

private:

	void pushBasePairToTokenizer(std::bitset<N> & tokenizer, char val) {
		tokenizer = ((tokenizer << bitsPerBasePair) & mask) | val;// getValueFromCharacter(val);
	}

	size_t bitsNeeded;
	size_t kmerWidth;
	std::bitset<N> mask;
	std::bitset<N> firstBasePairMask;

};

//This function is a "manual" implementation of the std::bitset Merizer.
//It seems that std::bitset comes with some overhead when try to do
//bit manipulation, so this class attempts to work around it by
//allow direct access to the bits used to encode kmers
template <typename T, size_t N> class Merizer<std::array<T, N>> {

	static_assert(std::is_unsigned<T>::value == true, "Array type must be unsigned");

public:

	static const size_t  bitsPerBasePair;// = 3;
	size_t const bitsPerSizeT = sizeof(T) * 8;
	static size_t neededBitsForKmerWidth(size_t kmerWidth) {
		return kmerWidth * bitsPerBasePair;
	};

	Merizer(size_t _kmerWidth) :kmerWidth(_kmerWidth) {

		bitsNeeded = neededBitsForKmerWidth(kmerWidth);

		if (bitsNeeded > N*bitsPerSizeT)
			throw KmerError(1, "Invalid bitset size for provided kmer width. Need at least " + std::to_string(bitsNeeded) + " bits");

		usedRegisters = (bitsNeeded - 1) / bitsPerSizeT + 1;

		fullBasePairsPerRegister = bitsPerSizeT / bitsPerBasePair;

	//For whatever reason, the below code gets compile incorrectly on linux release build. Windows release/debug and linux debug are fine.
	/*
		fullRegisterMask = 0;
		for (auto i = 0; i < bitsPerBasePair; i++) {
			fullRegisterMask |= (1 << (bitsPerSizeT - i - 1));
		}
	*/

		//This works.
        fullRegisterMask = (1 << bitsPerBasePair) - 1;
        fullRegisterMask <<= bitsPerSizeT - bitsPerBasePair;


		fullRegisterMaskShift = bitsPerSizeT - bitsPerBasePair;

		remainderBits = bitsNeeded - (usedRegisters - 1) * bitsPerSizeT;// bitsNeeded % bitsPerSizeT;
		lastRegisterMask = (1 << remainderBits) - 1;


	};

	//This ugliness was hidden by the stl::bitset implementation. Since multiple registers can
	//be used to make a "super register", it's necessary to make sure that overflows during bitshifting
	//are properly carried over to "high order" registers. 
	std::vector < std::array<T, N> > getMerTokens(std::vector<char> sequence) {

		std::vector<std::array<T, N>> rtn;

		if (sequence.size() < kmerWidth)
			return rtn;

		size_t totalKmers = sequence.size() - kmerWidth + 1;
		rtn.resize(totalKmers);

		std::array<T, N> tokenizer;
		for (auto & i : tokenizer)
			i = 0;

		//Need to prime tokenizer
		for (auto i = 0; i < kmerWidth - 1; i++)
			pushBasePairToTokenizer(tokenizer, sequence[i]);

		//Now start collecting tokens;
		for (auto i = 0; i < totalKmers; i++) {
			pushBasePairToTokenizer(tokenizer, sequence[i + kmerWidth - 1]);
			rtn[i] = tokenizer;
		}

		return rtn;

	};

	//Some more ugliness to decode an array into a sequence. Again, stl::bitset hid this nicely, but
	//doing it ourselves gives makes the difference between a slower and much faster implementation.
	std::vector<char> tokenToSequence(std::array<T, N> token) {

		std::vector<char> rtn(kmerWidth);

		auto bitsLeft = bitsPerSizeT;
		auto thisRegister = 0;
		for (auto i = 0; i < kmerWidth; i++) {

			auto currentElement = kmerWidth - 1 - i;

			auto readBits = std::min(bitsPerBasePair, bitsLeft);
			auto nextReadBits = bitsPerBasePair - readBits;

			T currentMask = (1 << readBits) - 1;
			T thisChar = token[thisRegister] & currentMask;
			token[thisRegister] >>= readBits;
			bitsLeft -= readBits;
			if (nextReadBits) {

				T nextMask = (1 << nextReadBits) - 1;
				thisRegister++;
				
				thisChar |= (token[thisRegister] & nextMask)<<readBits;
				token[thisRegister] >>= nextReadBits;
				bitsLeft = bitsPerSizeT - nextReadBits;
			}

			rtn[currentElement] = thisChar;// getCharacterFromValue(thisChar);
		}

		return rtn;

	};

private:

	void pushBasePairToTokenizer(std::array<T, N> & tokenizer, char val) {

		T incomingBits = val;// getValueFromCharacter(val);
		T highestBits = (tokenizer[0] & fullRegisterMask) >> fullRegisterMaskShift;
		tokenizer[0] = (tokenizer[0] << bitsPerBasePair) | incomingBits;

		for (auto i = 1; i < usedRegisters; i++) {
			incomingBits = highestBits;
			highestBits = (tokenizer[i] & fullRegisterMask) >> fullRegisterMaskShift;
			tokenizer[i] = (tokenizer[i] << bitsPerBasePair) | incomingBits;
		}

		tokenizer[usedRegisters - 1] &= lastRegisterMask;

	}

	size_t fullBasePairsPerRegister;
	size_t bitsNeeded;
	size_t kmerWidth;
	
	T lastRegisterMask;
	T fullRegisterMask;
	size_t fullRegisterMaskShift;
	size_t remainderBits;
	size_t usedRegisters;
	

};

//This is hanging out here due to some nuances in the C++ language when using static consts with templates.
template <typename T, size_t N> const size_t Merizer<std::array<T, N>>::bitsPerBasePair = 3;
