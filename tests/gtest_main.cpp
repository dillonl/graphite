#include "gtest/gtest.h"

#include "VCFFileTests.hpp"
#include "BamAlignmentReaderTests.hpp"
#include "RegionTests.hpp"
#include "FileTests.hpp"
#include "VCFFileTests.hpp"
#include "GSSWTests.hpp"
#include "GSSWGraphTests.hpp"
#include "AlleleTests.hpp"
#include "VariantsTest.hpp"
#include "CompoundVariantTests.hpp"
#include "FastaReferenceTests.hpp"

#include "GraphPathAlignmentTests.hpp"
#include "GraphPathRefTests.hpp"
#include "BamAlignmentTests.hpp"
#include "BamHeaderReaderTests.hpp"

GTEST_API_ int main(int argc, char** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
