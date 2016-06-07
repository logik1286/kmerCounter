#pragma once

#include "kmerError.h"
#include "fastqReader.h"
#include <iostream>
#include "kmerCounter.h"
#include <fstream>
#include <chrono>

namespace KmerUtils {

	KmerError getAllArguments(int argc, char * argv[], std::vector<std::string> & rtn);
	KmerError getArgument(int idx, int argc, char * argv[], std::string & rtn);
	int displayError(KmerError err, std::ostream & displayTo);
	void displayFastqEntry(FastqReader::FastqEntry & entry, std::ostream & displayTo);

	struct ProgramParams {

		std::string outputLog;
		std::ostream & outputTerminal;
		std::string filePath;
		size_t kmerWidth;
		size_t topCount;
		size_t threshold;

		ProgramParams(std::ostream & _outputTerminal);

	};

	//Main program. Specialized by counter
	template <class CounterType> KmerError countMers(ProgramParams & params) {

		try {

			std::ofstream outputStream;

			if (params.outputLog.size() != 0) {
				outputStream.open(params.outputLog.c_str(), std::ios_base::binary);
				if (outputStream.fail())
					return KmerError(1, "Error opening output log");
			};

			//Attempt to open fastqFile
			std::ifstream inputStream(params.filePath.c_str(), std::ios_base::binary); //keep this in scope, otherwise inputStream will be destroyed an will cause reader to have an invalid stream.
			if (!inputStream.is_open())
				return KmerError(1, "Error opening input file " + params.filePath);

			FastqReader reader(inputStream);

			FastqReader::FastqEntry entry;

			//Initialize kmer counter
			CounterType counter(params.kmerWidth);

			size_t count = 0;
			auto start = std::chrono::high_resolution_clock::now();

			KmerError err;

			//This loop reads in from a FASTQ file until its empty
			while (inputStream.good()) {

				err = reader.getEntry(entry);

				if (inputStream.eof())
					break;

				if (err.isError()) {
					return displayError(err, params.outputTerminal);
				}

				if (count % 10000 == 0) params.outputTerminal << "processed " << count << std::endl;

				err |= counter.addSequence(entry.sequence);
				if (err.isError()) return err;
				count++;
			}

			params.outputTerminal << "Getting mers" << std::endl;

			//Get kmers, apply threshold to filter out small values and improve sort performance
			typename CounterType::MerList outMer;

			err = counter.getTopMers(outMer, params.topCount, params.threshold);
			if (err.isError()) return err;

			auto end = std::chrono::high_resolution_clock::now();

			double duration = std::chrono::duration<double>(end - start).count();

			if (count != 0) {
				params.outputTerminal << "Total time [s] : " << duration << std::endl;
				params.outputTerminal << "Time per entry [ms] : " << duration / count * 1000.0 << std::endl;
			}
			else {
				return KmerError(1, "Nothing read.");
			}

			//Sort kmers by frequency. Partial_sort will only sort enough of the array to get the top N elements
			size_t topReportedCount = std::min(params.topCount, outMer.size());

			start = std::chrono::high_resolution_clock::now();
			std::partial_sort(outMer.begin(), outMer.begin() + topReportedCount, outMer.end(), CounterType::MerListSortByFrequency);
			end = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration<double>(end - start).count();
			//params.outputTerminal << "Sorting time [s] : " << duration << std::endl;

			//Complain if there weren't enough kmers to create a full top kmer list
			if (topReportedCount != params.topCount) {
				params.outputTerminal << "Warning. Specified top " << params.topCount << " of kmers, only " << topReportedCount << " were found." << std::endl;
			}

			auto saturatedCount = counter.getMaximumCount();
			bool foundSaturated = false;

			//Review top kmers. Check for saturated counters
			for (auto i = 0; i < topReportedCount; i++) {
				std::string seq;
				err = counter.decodeSequence(outMer[i].first, seq);
				if (err.isError()) return err;

				if (outMer[i].second == saturatedCount)
					foundSaturated = true;

				if (!outputStream.is_open())
					params.outputTerminal << seq << "," << +outMer[i].second << std::endl;
				else
					outputStream << seq << "," << +outMer[i].second << std::endl;
			}

			//Let user know they might want to rerun the application with larger accumulators
			if (foundSaturated)
				params.outputTerminal << "Warning, counters appear to be saturating. Consider increasing precision parameter." << std::endl;

		}
		catch (KmerError & e) {
			//Just in case something is thrown
			return e;
		}

		return KmerError();

	}

	//This function determines the type of counter to use in the program given the counter type and kmer width
	template <typename Precision> KmerError createProgram(char counterType, ProgramParams & params) {


		switch (counterType) {

		case 0: //Sort

			if (params.kmerWidth < 3)
				return countMers<SortMerizedKmerCounter<std::array<unsigned char, 1>, Precision>>(params);
			else if (params.kmerWidth < 6)
				return countMers<SortMerizedKmerCounter<std::array<unsigned short, 1>, Precision>>(params);
			else if (params.kmerWidth < 11)
				return countMers<SortMerizedKmerCounter<std::array<unsigned int, 1>, Precision>>(params);
			else if (params.kmerWidth < 22)
				return countMers<SortMerizedKmerCounter<std::array<unsigned int, 2>, Precision>>(params);
			else if (params.kmerWidth < 43)
				return countMers<SortMerizedKmerCounter<std::array<size_t, 2>, Precision>>(params);
			else
				return countMers<SortMerizedKmerCounter<std::vector<char>, Precision>>(params);
			break;

		case 1: //Map

			if (params.kmerWidth < 3)
				return countMers<HashMerizedKmerCounter<std::array<unsigned char, 1>, Precision>>(params);
			else if (params.kmerWidth < 6)
				return countMers<HashMerizedKmerCounter<std::array<unsigned short, 1>, Precision>>(params);
			else if (params.kmerWidth < 11)
				return countMers<HashMerizedKmerCounter<std::array<unsigned int, 1>, Precision>>(params);
			else if (params.kmerWidth < 22)
				return countMers<HashMerizedKmerCounter<std::array<unsigned int, 2>, Precision>>(params);
			else if (params.kmerWidth < 43)
				return countMers<HashMerizedKmerCounter<std::array<size_t, 2>, Precision>>(params);
			else
				return countMers<HashMerizedKmerCounter<std::vector<char>, Precision>>(params);

			break;

		case 2: //Unordered map

			if (params.kmerWidth < 3)
				return countMers<UnorderedHashMerizedKmerCounter<std::array<unsigned char, 1>, Precision, ArrayHasher<unsigned char, 1>>>(params);
			else if (params.kmerWidth < 6)
				return countMers<UnorderedHashMerizedKmerCounter<std::array<unsigned short, 1>, Precision, ArrayHasher<unsigned short, 1>>>(params);
			else if (params.kmerWidth < 11)
				return countMers<UnorderedHashMerizedKmerCounter<std::array<unsigned int, 1>, Precision, ArrayHasher<unsigned int, 1>>>(params);
			else if (params.kmerWidth < 22)
				return countMers<UnorderedHashMerizedKmerCounter<std::array<unsigned int, 2>, Precision, ArrayHasher<unsigned int, 2>>>(params);
			else if (params.kmerWidth < 43)
				return countMers<UnorderedHashMerizedKmerCounter<std::array<size_t, 2>, Precision, ArrayHasher<size_t, 2>>>(params);
			else
				return countMers<UnorderedHashMerizedKmerCounter<std::vector<char>, Precision, VectorHasher<char>>>(params);

			break;

		}


		return KmerError(1, "Invalid counter type");

	}

};
