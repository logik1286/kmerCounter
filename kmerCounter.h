#pragma once

#include "kmerError.h"
#include "merizer.h"
#include "keyCache.h"
#include <memory>
#include <map>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <locale>
#include <limits> 

//KmerCounter implements hashing/sorting algorithms independent from how kmers are expressed.
//With one exception - it accepts the original GTACN characters as a string and then converts it into
//01234. This could be done at the beginning of the Merizer, however I chose to have it done by the
//KmerCounter to force the Merizers to work with a simpler representation of GTACN characters and
//minimize duplicate code.
template <typename C = size_t> class KmerCounter {

	static_assert(std::is_unsigned<C>::value == true, "KmerCounter precision should be unsigned.");


public:

	C getMaximumCount(){
		return std::numeric_limits<C>::max();
	};

	typedef std::pair<std::vector<char>, C> MerListEntry;
	typedef std::vector<MerListEntry> MerList;

	static bool MerListSortByKey(MerListEntry & a, MerListEntry & b) {
		return a.first > b.first;
	}

	static bool MerListSortByFrequency(MerListEntry & a, MerListEntry & b) {
		return a.second > b.second;
	}

	KmerCounter(size_t _kmerWidth) :kmerWidth(_kmerWidth) {};
	
	virtual KmerError addSequence(std::string sequence) = 0;
	virtual KmerError getTopMers(MerList & mers, size_t num, C threshold) = 0;

	KmerError decodeSequence(std::vector<char> encoded, std::string & decoded) {

		decoded.resize(encoded.size());
		for (auto i = 0; i < encoded.size(); i++) {

			switch (encoded[i]) {
			case 0:
				decoded[i] = 'G';
				break;
			case 1:
				decoded[i] = 'T';
				break;
			case 2:
				decoded[i] = 'C';
				break;
			case 3:
				decoded[i] = 'A';
				break;
			case 4:
				decoded[i] = 'N';
				break;
			default:
				return KmerError(1, "Got invalid encoded value");
			}


		}

		return KmerError();

	}

	KmerError encodeSequence(std::string decoded, std::vector<char> & encoded) {

		//Force string to be upper case. Supposedly (from what I read) lower case characters are converted to upper case in this kind of application.

		std::transform(decoded.begin(), decoded.end(), decoded.begin(), ::toupper);

		encoded.resize(decoded.size());

		for (auto i = 0; i < encoded.size(); i++) {

			switch (decoded[i]) {
			case 'G':
				encoded[i] = 0;
				break;
			case 'T':
				encoded[i] = 1;
				break;
			case 'C':
				encoded[i] = 2;
				break;
			case 'A':
				encoded[i] = 3;
				break;
			case 'N':
				encoded[i] = 4;
				break;
			default:
				return KmerError(1, "Got invalid sequence character");
			}

		}

		return KmerError();

	};


	virtual ~KmerCounter() {}; //Need to declare virtual destructor for proper cleanup since pure virtual functions are in use

protected:

	size_t kmerWidth;

};

template <class T, typename C = size_t> class MerizierKmerCounter : public KmerCounter<C> {

public:

	MerizierKmerCounter(size_t kmerWidth, size_t _cacheOnEntries) : KmerCounter<C>(kmerWidth), merizer(kmerWidth), cache("kmerCache",true), cacheOnEntries(_cacheOnEntries) {};
	virtual KmerError addSequence(std::string sequence) = 0;
	virtual KmerError getTopMers(typename KmerCounter<C>::MerList & mers, size_t num, C threshold) = 0;
	
	virtual ~MerizierKmerCounter() {}; //need virtual destructor for proper cleanup. don't really need to do anything for this base class destructor

protected:

	virtual KmerError flushToCache() = 0;

	size_t cacheOnEntries;
	Merizer<T> merizer;
	KeyCache<T, C> cache;
};


template <class T, typename C = size_t, class L = std::greater<T>> class HashMerizedKmerCounter : public MerizierKmerCounter<T, C> {
	

public:

	HashMerizedKmerCounter(size_t kmerWidth, size_t cacheOnEntries = 10000000) : MerizierKmerCounter<T, C>(kmerWidth, cacheOnEntries) {};

	KmerError addSequence(std::string sequence) {

		std::vector<char> encodedSequence;
		auto err = this->encodeSequence(sequence, encodedSequence);
		if (err.isError()) return err;
		auto tokens = this->merizer.getMerTokens(encodedSequence);

		for (auto & i : tokens) {
			auto & val = hashTable[i];
			if(val != std::numeric_limits<C>::max())
				val++;
		};

		if(hashTable.size() > this->cacheOnEntries)
			return flushToCache();

		return KmerError();

	};

	KmerError getTopMers(typename KmerCounter<C>::MerList & mers, size_t num, C threshold = 1) {


        KmerError err = flushToCache();
        if(err.isError()) return err;

        std::vector<typename KeyCache<T,C>::EntryType> topEntries;

        err = this->cache.getTopKeys(topEntries, num, threshold);
		if(err.isError()) return err;
		mers.clear();
        mers.reserve(topEntries.size());

        for(auto i = 0; i < topEntries.size(); i++){
            auto const & seq = this->merizer.tokenToSequence(topEntries[i].first);
            mers.push_back(std::make_pair(seq, topEntries[i].second));
        };

		return KmerError();

	};

protected:

	KmerError flushToCache(){

		KmerError err;
       typename KeyCache<T,C>::EntryType entry;

		for(auto & i : hashTable){
			entry.first = i.first;
			entry.second = i.second;

			err = this->cache.incrementKey(entry);
			if(err.isError()) return err;
		};

		hashTable.clear();


		return KmerError();
	};

private:

	std::map<T, C, L> hashTable;

};

template <typename T, size_t N> struct ArrayHasher
{
	T operator()(const std::array<T, N> & k) const
	{
		auto temp = std::hash<T>()(k[0]);

		for (auto i = 1; i < N; i++)
			temp ^= std::hash<T>()(k[i]);

		return temp;

	}
};

//where is this from?
template <typename T> struct VectorHasher
{
	std::size_t operator()(const std::vector<T> & vec) const
	{
		std::size_t seed = vec.size();
		for (auto& i : vec) {
			seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		return seed;

	}
};



template <class T, typename C = size_t, class L = std::greater<T>> class SortMerizedKmerCounter : public MerizierKmerCounter<T, C> {


	public:

		SortMerizedKmerCounter(size_t kmerWidth, size_t cacheOn = 100000) : MerizierKmerCounter<T, C>(kmerWidth, cacheOn) {
			tokens.reserve(this->cacheOnEntries);
		};

		KmerError addSequence(std::string sequence) {

			std::vector<char> encodedSequence;
			auto err = this->encodeSequence(sequence, encodedSequence);
			if (err.isError()) return err;
			tokens.push_back(std::make_shared<std::vector<T>>(this->merizer.getMerTokens(encodedSequence)));

			if(tokens.size() > this->cacheOnEntries)
				return flushToCache();

			return KmerError();

		};

        KmerError getTopMers(typename KmerCounter<C>::MerList & mers, size_t num, C threshold = 1) {

            KmerError err = flushToCache();
            if(err.isError()) return err;

            std::vector<typename KeyCache<T,C>::EntryType> topEntries;

            err = this->cache.getTopKeys(topEntries, num, threshold);

            mers.reserve(topEntries.size());

            for(auto i = 0; i < topEntries.size(); i++){
                auto const & seq = this->merizer.tokenToSequence(topEntries[i].first);
                mers.push_back(std::make_pair(seq, topEntries[i].second));
            };


            return KmerError();

        };

#if 0
		KmerError getTopMers(typename KmerCounter<C>::MerList & mers, size_t num, C threshold = 1) {

			mers.clear(); //clean everything out

			KmerError err = flushToCache();
			if(err.isError()) return err;

			//this->cache.dumpAllKeys();

			std::vector<typename KeyCache<T,C>::EntryType> topEntries;

			err = this->cache.getTopKeys(topEntries, num, threshold);
			if(err.isError()) return err;

			mers.reserve(topEntries.size());

			for(auto i = 0; i < topEntries.size(); i++){
				auto const & seq = this->merizer.tokenToSequence(topEntries[i].first);
				mers.push_back(std::make_pair(seq, topEntries[i].second));
			};

			return KmerError();

		};
#endif
	protected:

		KmerError flushToCache(){

			typename KeyCache<T,C>::EntryType entry;

			//copy all tokens into single vector
			size_t totalTokens = 0;
			for (auto i = 0; i < tokens.size(); i++) {
				totalTokens += tokens[i]->size();
			}

			if (totalTokens == 0)
				return KmerError(1, "No mers generated.");

			std::vector<T> tempArray; tempArray.reserve(totalTokens);
			size_t copyCount = 0;
			for (auto i = 0; i < tokens.size(); i++) {
				tempArray.insert(tempArray.begin() + copyCount, tokens[i]->begin(), tokens[i]->end());
				copyCount += tokens[i]->size();
			}

			//Sort tokens, consolidate
			std::sort(tempArray.begin(), tempArray.end(), L());

			T thisToken = tempArray[0];
			C count = 0;

			KmerError err;

			for (auto i = 0; i < tempArray.size(); i++) {

				if (tempArray[i] != thisToken) {
					//create new entry if above threshold
					entry.first = thisToken;
					entry.second = count;
					err = this->cache.incrementKey(entry);
					if(err.isError()) return err;
					thisToken = tempArray[i];
					count = 0;
				}

				if (count != std::numeric_limits<C>::max())
					count++;

			}

			if(count!=0){
				//write out the remainder
				entry.first = thisToken;
				entry.second = count;

				err = this->cache.incrementKey(entry);
				if(err.isError()) return err;


			}


			tokens.clear();

			return KmerError();
		};


	private:

		size_t expectedEntries;
		float resizeFactor;
		std::vector < std::shared_ptr<std::vector<T>>> tokens;

};


template <class T, typename C = size_t, class L = std::hash<T>> class UnorderedHashMerizedKmerCounter : public MerizierKmerCounter<T, C> {
	//add static assert that counter precision supports
	public:

		UnorderedHashMerizedKmerCounter(size_t kmerWidth, size_t cacheOn = 100000000) : MerizierKmerCounter<T, C>(kmerWidth, cacheOn) {};

		KmerError addSequence(std::string sequence) {

			std::vector<char> encodedSequence;
			auto err = this->encodeSequence(sequence, encodedSequence);
			if (err.isError()) return err;
			auto tokens = this->merizer.getMerTokens(encodedSequence);

			for (auto & i : tokens) {
				auto & val = hashTable[i];
				if(val != std::numeric_limits<C>::max())
					hashTable[i]++;
			};

			if(hashTable.size()>this->cacheOnEntries)
				return flushToCache();

			return KmerError();

		};

		KmerError getTopMers(typename KmerCounter<C>::MerList & mers, size_t num, C threshold = 1) {

			KmerError err = flushToCache();
			if(err.isError()) return err;

			std::vector<typename KeyCache<T,C>::EntryType> topEntries;

			err = this->cache.getTopKeys(topEntries, num, threshold);

			mers.reserve(topEntries.size());

			for(auto i = 0; i < topEntries.size(); i++){
				auto const & seq = this->merizer.tokenToSequence(topEntries[i].first);
				mers.push_back(std::make_pair(seq, topEntries[i].second));
			};


			return KmerError();

		};

	protected:

		KmerError flushToCache(){

			KmerError err;

			typename KeyCache<T,C>::EntryType entry;

			for(auto & i : hashTable){
				entry.first = i.first;
				entry.second = i.second;

				err = this->cache.incrementKey(entry);
				if(err.isError()) return err;
			};

			hashTable.clear();


			return KmerError();
		};


	private:


		std::unordered_map<T, C, L> hashTable;

};
