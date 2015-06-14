#include <stdio.h>

#include "gtest/gtest.h"

/*
#include "BamAlignmentTests.hpp"
#include "VCFFileTests.hpp"
#include "ReferenceTest.hpp"
#include "VariantListVCFPreloadedTests.hpp"
#include "plugins/TestIncludes.h"
#include "FileTests.hpp"
*/

/*
#include "SequenceManagerTests.hpp"
#include "BuildGraphTests.hpp"
*/
#include "VariantsTest.hpp"
#include "RegionTests.hpp"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
