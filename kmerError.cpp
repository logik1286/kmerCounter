#include "kmerError.h"

using namespace std;

KmerError::KmerError() {

}

KmerError::KmerError(const std::exception & err) {
	
	errs.push_back(make_pair(defaultError, err.what()));

}


KmerError::KmerError(size_t idx) {

	errs.push_back(make_pair(idx, "no error"));

}

KmerError::KmerError(size_t idx, std::string what) {

	errs.push_back(make_pair(idx, what));

}

bool KmerError::isError() {

	return errs.size() != 0;

}

KmerError & KmerError::operator |= (const KmerError & rhs) {
	errs.insert(errs.end(), rhs.errs.begin(), rhs.errs.end());
	return *this;
}

KmerError & KmerError::operator |= (const exception & rhs) {

	errs.push_back(make_pair(defaultError, rhs.what()));
	return *this;
}

KmerError & KmerError::operator = (const KmerError & rhs) {

	errs = rhs.errs;
	return *this;
}


KmerError & KmerError::operator = (const exception & rhs) {

	errs.clear();
	errs.insert(errs.end(), make_pair(defaultError, rhs.what()));
	return *this;
}


KmerError::ErrorList KmerError::getErrorDescriptions() {
	return errs;
}

KmerError::ErrorEntry KmerError::getLastErrorDescription() {
	
	if (errs.size() == 0) {
		return ErrorEntry(); //blank entry
	}

	return *errs.rbegin();


};

const char* KmerError::what() const throw () {

	if (errs.size())
		return "Program errors detected. Retrieve with getErrorDescriptions";

	return "no error";

};

KmerError::~KmerError() throw(){};
