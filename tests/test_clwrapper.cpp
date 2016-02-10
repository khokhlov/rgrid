#define CATCH_CONFIG_MAIN


#include "catch.hpp"
#include "rgrid/clwrapper.h"

#include <iostream>

using namespace std;

TEST_CASE(
		"CLWrapper",
		"constructor"
	 )
{
	clwrapper::CLWrapper::instance();
}
