#ifndef GWIZ_TESTS_BAMALIGNMENTTESTS_HPP
#define GWIZ_TESTS_BAMALIGNMENTTESTS_HPP

#include "core/alignments/BamAlignmentReader.h"
#include "TestConfig.h"

#include "api/BamReader.h"
#include "api/BamAlignment.h"

TEST(BamAlignmentTests, OpenBamAlignmentInvalidTest)
{
    bool success = false;
	try
	{
		std::string path = TEST_INVALID_FILE;
		auto reader = std::make_shared< gwiz::BamAlignmentReader >(path);
	}
	catch (...)
	{
        success = true;
	}
	EXPECT_TRUE(success);
}

TEST(BamAlignmentTests, OpenBamAlignmentValidTest)
{
    bool success = true;
	try
	{
		std::string path = TEST_BAM_FILE;
		auto reader = std::make_shared< gwiz::BamAlignmentReader >(path);
	}
	catch (...)
	{
        success = false;
	}
	EXPECT_TRUE(success);
}

TEST(BamAlignmentTests, ReadOneBamAlignmentTest)
{
	std::string path = TEST_BAM_FILE;
	auto reader = std::make_shared< gwiz::BamAlignmentReader >(path);
	// EXPECT_TRUE(success);
}

#endif //GWIZ_TESTS_BAMALIGNMENTTESTS_HPP
