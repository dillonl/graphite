#include "gtest/gtest.h"

#include "FileTests.hpp"
#include "VCFFileTests.hpp"
#include "SequenceManagerTests.hpp"
#include "VariantsTest.hpp"
#include "RegionTests.hpp"
#include "BamAlignmentReaderTests.hpp"
#include "BamAlignmentManagerTests.hpp"
#include "GSSWTests.hpp"
#include "plugins/TestIncludes.h"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
