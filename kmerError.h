#pragma once

#include <exception>
#include <string>
#include <vector>

//Error class. So useful.
class KmerError : public std::exception {

	const size_t noError = 0;
	const size_t defaultError = 1;
	typedef size_t ErrorCode;
	typedef std::pair<ErrorCode, std::string> ErrorEntry;
	typedef std::vector<ErrorEntry> ErrorList;

public:

	KmerError();
	KmerError(const std::exception & err);
	KmerError(size_t idx);
	KmerError(size_t idx, std::string what);

	bool isError();

	KmerError & operator |= (const KmerError & rhs);
	KmerError & operator |= (const std::exception & rhs);
	KmerError & operator = (const std::exception & rhs);
	KmerError & operator = (const KmerError & rhs);

	ErrorList getErrorDescriptions();
	ErrorEntry getLastErrorDescription();

	const char* what() const throw ();

	~KmerError() throw();

private:

	ErrorList errs;

};
