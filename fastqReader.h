#pragma once

#include "kmerError.h"
#include <istream>

//Class that manages reading from Fastq file
class FastqReader {

public:

	FastqReader(std::istream & inputStream);

	struct FastqEntry {
		std::string seqID;
		std::string sequence;
		std::string alignment;
		std::string quality;
	};

	KmerError getEntry(FastqEntry & entry);

private:

	KmerError checkSequence(std::string & checkSequence);

	std::istream & inputStream;
	
};