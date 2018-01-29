/*
 * Test.cpp
 *
 *  Created on: Jun 8, 2014
 *      Author: alinsch
 *
 *	A small program to test the error handling object
 */
#include <iostream>
#include "gslpp/error_handling/Error.h"
#include "gslpp/error_handling/test/Test.h"

namespace gslpp{
namespace error_handling{

void RunTest::run_test(){
	std::cout << "\n\nTest of the error handling methods" <<std::endl;

	std::cout << "Test success" <<std::endl;
}

};
};
