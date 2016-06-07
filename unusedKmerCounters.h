#pragma once

#include "kmerCounter.h"
#if 0
template <class T, typename C = size_t, class L = std::greater<T>> class FlyingSortMerizedKmerCounter : public MerizierKmerCounter<T, C> {
	//add static assert that counter precision supports
public:

	FlyingSortMerizedKmerCounter(size_t kmerWidth, size_t _expectedEntries = 1000000) : MerizierKmerCounter<T>(kmerWidth), entriesBeforeHash(_expectedEntries) {
		tokens.reserve(entriesBeforeHash);
	};

	KmerError addSequence(std::string sequence) {

		std::vector<char> encodedSequence;
		auto err = this->encodeSequence(sequence, encodedSequence);
		if (err.isError()) return err;
		tokens.push_back(this->merizer.getMerTokens(encodedSequence));

		if (tokens.size() == entriesBeforeHash)
			return consolidateTokens();

		return KmerError();

	};

	KmerError consolidateTokens() {

		size_t totalTokens = 0;

		for (auto & v : tokens)
			totalTokens += v.size();

		if (totalTokens == 0) return KmerError();

		std::vector<T> temp; temp.reserve(totalTokens);

		for (auto & v : tokens) {
			temp.insert(temp.end(), v.begin(), v.end());
		};

		std::vector<std::pair<T, size_t>> tempPair; tempPair.resize(temp.size());

		for (auto i = 0; i < totalTokens; i++) {
			tempPair[i].first = temp[i];
			tempPair[i].second = 1;
		}

		auto tempCount = countedTokens;


		tempCount.insert(tempCount.end(), tempPair.begin(), tempPair.end());
		std::sort(tempCount.begin(), tempCount.end());

		countedTokens.clear();

		T thisToken = tempCount[0].first;
		size_t count = 0;
		for (auto i = 0; i < tempCount.size(); i++) {

			if (tempCount[i].first != thisToken) {
				countedTokens.push_back(std::make_pair(thisToken, count));
				thisToken = tempCount[i].first;
				count = 0;
			}

			count += tempCount[i].second;

		}


		tokens.clear();

		return KmerError();
	};

	KmerError getMers(MerList & mers, size_t threshold = 1) {

		mers.clear(); //clean everything out

		auto err = consolidateTokens();

		if (err.isError()) return err;

		mers.resize(countedTokens.size());

		for (auto i = 0; i < mers.size(); i++) {

			auto & token = countedTokens[i];
			if (token.second > threshold) {
				auto & seq = this->merizer.tokenToSequence(token.first);
				mers[i] = std::make_pair(seq, token.second);
			}
		}

		return KmerError();

	};

private:

	size_t entriesBeforeHash;
	std::vector <std::vector<T>> tokens;
	std::vector<std::pair<T, size_t>> countedTokens;

};



template <typename C> class TreeKmerCounter : public KmerCounter<C> {

public:


	TreeKmerCounter(size_t _kmerWidth) : KmerCounter(_kmerWidth) {


	}

	KmerError addSequence(std::string sequence) {

		//get all mers
		std::vector<char> encodedSequence;
		auto err = encodeSequence(sequence, encodedSequence);
		if (err.isError()) return err;

		std::vector<char> working(kmerWidth);

		for (auto i = encodedSequence.begin(); i + kmerWidth != encodedSequence.end(); ++i) {
			working.assign(i, i + kmerWidth);
			err |= descendTree(working, working.begin(), root);
		}

		return err;

	}

	KmerError getMers(MerList & mers, size_t threshold = 1) {

		size_t totalMers = root.accessCount;
		mers.reserve(totalMers);
		vector<char> working(kmerWidth); //allocate temporary buffer

		return parseTree(mers, working, 0, root, threshold);

	}


private:

	class Node {

	public:

		Node::Node() {
			accessCount = 0;
		}

		Node & Node::operator [](char & val) {

			accessCount++;
			return children[val];

		}

		std::map<char, Node> children;
		C accessCount;

	};

	Node root;


	KmerError  descendTree(std::vector<char>  & seq, std::vector<char>::iterator it, Node & current) {

		if (it == seq.end()) {
			current.accessCount++;
			return KmerError();
		}

		char thisBP = *it;
		Node & next = current[thisBP];

		return descendTree(seq, ++it, next);


	}

	KmerError  parseTree(MerList & workingList, std::vector<char> & workingSeq, size_t pos, Node & current, C threshold) {

		if (pos == workingSeq.size()) {
			//creaet new entry
			if (current.accessCount != 0 && current.accessCount > threshold)
				workingList.push_back(make_pair(workingSeq, current.accessCount));
			return KmerError();
		}

		KmerError err;

		for (auto & i : current.children) {
			workingSeq[pos] = i.first;
			err |= parseTree(workingList, workingSeq, pos + 1, i.second);
		};

		return err;

	};


};
#endif
