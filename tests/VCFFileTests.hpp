#ifndef GWIZ_VCF_FILE_TESTS_HPP
#define GWIZ_VCF_FILE_TESTS_HPP

#include "gtest/gtest.h"

#include <stdexcept>

#include "core/variants/VCFFileReader.h"
#include "core/Region.h"

#include "TestConfig.h"

TEST(VCFFileReaderTests, TestVCFFileReaderGetRegion)
{
	std::string chrom = "Y";
	gwiz::position startPosition = 2851690;
	gwiz::position endPosition = 2900741;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	gwiz::Region::SharedPtr region = std::make_shared< gwiz::Region >(regionString);
	std::string path = TEST_VCF_FILE;
	auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(path);
	vcfFileReader->Open(region);
	gwiz::Variant::SharedPtr variantPtr;
	vcfFileReader->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), startPosition);
	while (vcfFileReader->getNextVariant(variantPtr))
	{
		ASSERT_STREQ(chrom.c_str(), variantPtr->getChrom().c_str());
		ASSERT_LE(startPosition, variantPtr->getPosition());
		ASSERT_GE(endPosition, variantPtr->getPosition());
	}
	ASSERT_EQ(variantPtr->getPosition(), endPosition);
}

TEST(VCFFileReaderTests, TestVCFFileReaderGetRegionBeforeStart)
{
	std::string chrom = "Y";
	gwiz::position startPosition = 2655179;
	gwiz::position endPosition = 2656000;
	gwiz::position firstPosition = 2655180;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	gwiz::Region::SharedPtr region = std::make_shared< gwiz::Region >(regionString);
	std::string path = TEST_VCF_FILE;
	auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(path);
	vcfFileReader->Open(region);
	gwiz::Variant::SharedPtr variantPtr;
	vcfFileReader->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
}

TEST(VCFFileReaderTests, TestVCFFileReaderGetRegionEnd)
{
	std::string chrom = "Y";
	gwiz::position startRegionPosition = 28760931;
	gwiz::position endRegionPosition = 38770931;
	gwiz::position endPosition = 28770931;
    std::string regionString = chrom + ":" + std::to_string(startRegionPosition) + "-" + std::to_string(endRegionPosition);
	gwiz::Region::SharedPtr region = std::make_shared< gwiz::Region >(regionString);
	std::string path = TEST_VCF_FILE;
	auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(path);
	vcfFileReader->Open(region);
	gwiz::Variant::SharedPtr variantPtr;
	while (vcfFileReader->getNextVariant(variantPtr))
	{
	}
	ASSERT_EQ(variantPtr->getPosition(), endPosition);
}

/*
 * Makes sure that a region without any variants throws an exception
 */
TEST(VCFFileReaderTests, TestVCFFileReaderGetNoRegion)
{
	std::string chrom = "Y";
	gwiz::position startRegionPosition = 0;
	gwiz::position endRegionPosition = 200;
    std::string regionString = chrom + ":" + std::to_string(startRegionPosition) + "-" + std::to_string(endRegionPosition);
	gwiz::Region::SharedPtr region = std::make_shared< gwiz::Region >(regionString);
	std::string path = TEST_VCF_FILE;
	auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(path);
	ASSERT_ANY_THROW(vcfFileReader->Open(region));
}

#endif // GWIZ_VCF_FILE_TESTS_HPP
