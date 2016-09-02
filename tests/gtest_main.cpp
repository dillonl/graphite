#include "gtest/gtest.h"

#include "RegionTests.hpp"
#include "VariantsTest.hpp"
/*
#include "FileTests.hpp"
#include "VCFFileTests.hpp"
#include "SequenceManagerTests.hpp"
#include "CompoundVariantTests.hpp"
#include "BamAlignmentReaderTests.hpp"
#include "BamAlignmentManagerTests.hpp"
#include "GSSWTests.hpp"
#include "AdjudicationTests.hpp"
#include "AlleleTests.hpp"
#include "plugins/TestIncludes.h"
*/
// #include "S3FileTests.hpp"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
