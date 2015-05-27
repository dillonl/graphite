#ifndef GWIZ_TESTS_BAMALIGNMENTTESTS_HPP
#define GWIZ_TESTS_BAMALIGNMENTTESTS_HPP

#include "core/alignment/BamAlignmentReader.h"
#include "core/alignment/BamAlignment.h"
#include "TestConfig.h"

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
	bool success = false;
	std::string path = TEST_BAM_FILE;
	auto reader = std::make_shared< gwiz::BamAlignmentReader >(path);

	gwiz::IAlignment::SharedPtr bamAlignmentPtr;
	if (reader->getNextAlignment(bamAlignmentPtr)) // not a great test but we'll add better one's later
	{
		success = bamAlignmentPtr.get() != NULL;
	}

	EXPECT_TRUE(success);
}

#endif //GWIZ_TESTS_BAMALIGNMENTTESTS_HPP
