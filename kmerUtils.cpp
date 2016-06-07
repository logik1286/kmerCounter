#include "kmerUtils.h"

using namespace std;

KmerUtils::ProgramParams::ProgramParams(std::ostream & _outputTerminal) :outputTerminal(_outputTerminal) {};

KmerError KmerUtils::getAllArguments(int argc, char * argv[], vector<string> & rtn) {

	if (argc <= 0)
		return KmerError(1, "number of command arguments must be greater than or equal to zero");

	rtn.resize(argc);

	KmerError err;

	for (auto i = 0; i < argc; i++) {

		err |= getArgument(i, argc, argv, rtn[i]);

	}

	return err;

}

KmerError KmerUtils::getArgument(int idx, int argc, char * argv[], string & rtn) {

	if (idx > argc)
		return KmerError(1, "invalid argument index. provided" + to_string(idx) + ", maximum is " + to_string(argc - 1));

	try {
		rtn = string(argv[idx]);
	}
	catch (std::exception & err) {
		return err;
	}

	return KmerError();

}

int KmerUtils::displayError(KmerError err, std::ostream & displayTo){

	auto desc = err.getErrorDescriptions();

	if (desc.size() == 0)
		return 0;

	displayTo << "Got errors:" << endl;

	for (auto & i : desc) {
		displayTo << i.first << " : " << i.second << endl;
	}

	return 1;

}

void KmerUtils::displayFastqEntry(FastqReader::FastqEntry & entry, std::ostream & displayTo) {

	displayTo << "----FASTQ Entry----" << endl;
	displayTo << "Sequence ID: " << entry.seqID << endl;
	displayTo << "Sequence :" << entry.sequence << endl;
	displayTo << "Aligment :" << entry.alignment << endl;
	displayTo << "Qualiy : " << entry.quality << endl;

	return;

}
