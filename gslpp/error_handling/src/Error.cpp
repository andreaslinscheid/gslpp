/*
 * Error.cpp
 *
 *  Created on: May 31, 2014
 *      Author: alinsch
 */
#include <string>
#include <cstdlib>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include "gslpp/error_handling/Error.h"

namespace gslpp {
namespace error_handling {

void Error::call_error_with_code_ref(int line,
		std::string const& file,
		std::string const& description,
		int errorCode){
	std::cout << "ERROR: Occurred in file " << file << " line " << line <<": " << std::endl;
	std::cout << "\tProblem description : " <<description << std::endl;

	//If we are using debugging symbols generate a stack trace
#ifndef __NDEBUG
	handler(errorCode);
#endif
	std::exit(errorCode);
}
void Error::call_error_no_ref(std::string const& description,int errorCode){
	std::cout << "ERROR occurred ";
	std::cout << "\tProblem description : " <<description << std::endl;

	//If we are using debugging symbols generate a stack trace
#ifdef __DEBUG
	handler(errorCode);
#endif
	std::exit(errorCode);
}

//backtrace and backtrace_symbols_fd do not malloc and should thus be ok to call on throw.
void signal_handler(int signal) {
	 void * stackEntries[20];
	  size_t size;
	  //
	  //get pointers to all entries on the stack
	  size = backtrace(stackEntries, 20);
	  //
	  // print out all the frames to stderr
	  std::cerr << "Error: signal %d:\n" << signal << std::endl;
	  backtrace_symbols_fd(stackEntries, size, STDERR_FILENO);
	  std::exit(1);
}

void Error::handler(int signal) {
	signal_handler(signal);
}

void Error::generate_stacktrace_on_segfault() {
	//
	//install signal handler
	signal(SIGSEGV, signal_handler);
}

}; /* namespace error_handling */
}; /* namespace gslpp */

