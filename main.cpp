#include "fastqReader.h"
#include "kmerError.h"
#include <fstream>

#include "kmerUtils.h"
#include "kmerCounter.h"
#include <algorithm>
#include <chrono>

using namespace KmerUtils;
using namespace std;

//Main entry point
int main(int argc, char * argv[]) {

	ProgramParams params(cout);

	vector<string> args;
	KmerError err = getAllArguments(argc, argv, args);
	if (err.isError()) return displayError(err, params.outputTerminal);


	//Check argument size
	if (args.size() != 6 && args.size() != 7) {

		params.outputTerminal << "Invalid number of arguments. Usage:" << endl;
		params.outputTerminal << "kmerCounter <inputFile> <kmerSize> <topKmersToReport> <precision> <counterType> [output]" << endl;
		params.outputTerminal << "inputFile : FASTQ file to process" << endl;
		params.outputTerminal << "kmerSize  : Number of base pairs in a mer" << endl;
		params.outputTerminal << "topKmersToReport : The number of most frequent kmers to report" << endl;
		params.outputTerminal << "precision : size of accumulators. 0 = 1 byte (max 2^8-1), 1 = 2 bytes (max 2^16-1), 2 = 4 bytes (max 2^32-1), 3 = 8 bytes (max 2^64-1)" << endl;
		params.outputTerminal << "counterType: sorting algorithm to use. 0 = Sort and accumulate (fast, worse memory), 1 = Ordered hash Map (slow, good memory, consistent performance), 2 = Unordered hash map (better speed, good memory, delays during rehashing)" << endl;
		params.outputTerminal << "output [optional] : output file to write top kmers to" << endl;

		params.outputTerminal << endl;

		return 1;

	};

	char precision = 0;
	char counterType = 0;

	//Get arguments
	try {

		params.filePath = args[1];
		params.kmerWidth = stoul(args[2]);

		if(params.kmerWidth == 0)
			return displayError(KmerError(1, "Invalid kmer width :" + to_string(params.kmerWidth)), params.outputTerminal);


		params.topCount = stoul(args[3]);

        if(params.topCount == 0)
            return displayError(KmerError(1, "Invalid topKmersToReport:" + to_string(params.topCount)), params.outputTerminal);

		precision = stoul(args[4]);
		counterType = stoul(args[5]);
		params.threshold = 0; //deprecated
		if (args.size() == 7)
			params.outputLog = args[6];

	}
	catch (std::exception & e) {
		return displayError(e, params.outputTerminal);
	}

	//specialize program per counter precision. createProgram defined in kmerUtils.h
	switch (precision) {
	case 0:
		err = createProgram<unsigned char>(counterType, params);
		break;
	case 1:
		err = createProgram<unsigned short>(counterType, params);
		break;
	case 2:
		err = createProgram<unsigned int>(counterType, params);
		break;
	case 3:
		err = createProgram<size_t>(counterType, params);
		break;
	default:
		return displayError(KmerError(1, "Invalid precision"), params.outputTerminal);
	}


	return displayError(err, params.outputTerminal);

};
