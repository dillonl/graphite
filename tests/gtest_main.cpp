#include "gtest/gtest.h"

#include "VCFFileTests.hpp"
#include "BamAlignmentReaderTests.hpp"
#include "RegionTests.hpp"
#include "FileTests.hpp"
#include "VCFFileTests.hpp"
#include "SequenceManagerTests.hpp"
#include "GSSWTests.hpp"
#include "AlleleTests.hpp"


/*
#include "VariantsTest.hpp" // segfault

#include "BamAlignmentManagerTests.hpp"
#include "CompoundVariantTests.hpp"
#include "AdjudicationTests.hpp"


#include "plugins/TestIncludes.h"
*/
// #include "S3FileTests.hpp"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
