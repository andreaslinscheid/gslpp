/*
 * Error.h
 *
 *  Created on: May 31, 2014
 *      Author: Andreas Linscheid
 */

#ifndef GSLPP_ERROR_HANDLING_ERROR_H_
#define GSLPP_ERROR_HANDLING_ERROR_H_

#include <string>
#include <typeinfo>

namespace gslpp {
namespace error_handling {

/**
 * Quit the execution with an error message and possibly debug info.
 * If __NDEBUG is not defined this will cause to generate a stacktrace.
 */
class Error {
public:
	typedef enum {
		FILE_IO_ERR = 1,
		OUT_OF_BOUNDS = 10,
		ACCESS_WITHOUT_INIT = 11,
		INTERNAL_LOGIC_CHECK_FAILED = 43,
		INPUT_ERROR = 100
	}ERR_TYPE;

	//ctor that is not calling an error
	Error(){};

	//install a signal handler on a segfault causing a stactrace to be generated
	void generate_stacktrace_on_segfault();

	//ctor that is calling an error directly
	Error(int line,std::string const& file,std::string const& description,int errorCode){
		call_error_with_code_ref(line,file,description,errorCode);
	}

	//ctor that is calling an error directly without the reference to where it appeared
	Error(std::string const& description,int errorCode){
		call_error_no_ref(description,errorCode);
	}

	//ctor that is calling an error directly
	Error(int line,std::string const& file,std::string const& description,size_t errorCode){
		call_error_with_code_ref(line,file,description,static_cast<int>(errorCode));
	}
private:

	//generate a stack trace for debugging. Works with glibc, though may not be portable.
	void handler (int signal);

	//call error including the line and the file where the error was called from
	void call_error_with_code_ref(int line,std::string const& file,std::string const& description,int errorCode);

	//call error with no reference to where the error occurred.
	void call_error_no_ref(std::string const& description,int errorCode);
};

}; /* namespace error_handling */
}; /* namespace gslpp */
#endif /* GSLPP_ERROR_HANDLING_ERROR_H_ */
