#include "fastqReader.h"
#include <algorithm>
#include <locale>

using namespace std;

FastqReader::FastqReader(std::istream & _path):inputStream(_path) {
}

//Gets entry, and ensures that entry isn't corrupted. This does assume that there are no blank lines between entries or within a record.
KmerError FastqReader::getEntry(FastqEntry & entry) {

	//check input stream
	if (inputStream.fail())
		return KmerError(5, "Input stream is bad. Is it open?");

	//Reset entry
	entry = FastqEntry();

	//Read entry info
	getline(inputStream, entry.seqID);
	getline(inputStream, entry.sequence);
	getline(inputStream, entry.alignment);
	getline(inputStream, entry.quality);

	if (inputStream.fail())
		return KmerError(1, "Error while reading from FASTQ file");
	
	try {

		auto err = checkSequence(entry.sequence);
		if (err.isError()) return err;

		if (entry.seqID[0] != '@')
			return KmerError(2, "Malformed FASTQ entry. Sequence ID doesn't begin with @");

		if (entry.alignment[0] != '+')
			return KmerError(3, "Malformed FASTQ entry. Third line doesn't begin with +");

	}
	catch (std::exception & e) {
		return KmerError(4, "Error checking entry. Possible blank lines in FASTQ entry?" + string(e.what()));
	}

	return KmerError();

}

//Make sure that the sequence doesn't contain bogus characters.
KmerError FastqReader::checkSequence(std::string & sequence) {

	std::locale loc;

	for (auto & i : sequence) {

		char t = std::toupper(i, loc); //upper or lower case is valid, convert to upper and only check that

		if (t != 'G' && t != 'A' && t != 'T' && t != 'C' && t != 'N')
			return KmerError(1, "Got invalid base pair with ASCII value : " + to_string(i));

	};

	return KmerError();

}