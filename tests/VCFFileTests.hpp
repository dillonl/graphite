#ifndef GRAPHITE_VCF_FILE_TESTS_HPP
#define GRAPHITE_VCF_FILE_TESTS_HPP

#include <stdexcept>

#include "core/variant/IVariant.h"
#include "core/variant/VCFManager.h"
#include "core/region/Region.h"

#include "TestConfig.h"

TEST(VCFFileReaderTests, ReadAllChrom20)
{
	std::string chrom = "20";
    std::string regionString = chrom;
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	uint32_t totalCount = 181;
	std::string path = TEST_VCF_FILE;
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 0;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		ASSERT_STREQ(variantPtr->getChrom().c_str(), chrom.c_str());
		++count;
	}
	ASSERT_EQ(count, totalCount); // 181 is the number of variants in chrom 20 for the TEST_VCF_FILE
}

TEST(VCFFileReaderTests, ReadRegionAtStartNonExactBeginingAndEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 1;
	graphite::position endPosition = 4000000;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 3728155;
	uint32_t totalCount = 12;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 1;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
		++count;
	}
	ASSERT_EQ(tmpPosition, lastPosition);
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadRegionAtStartExactBeginingNonExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 909434;
	graphite::position endPosition = 4000000;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 3728155;
	uint32_t totalCount = 12;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 1;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
		++count;
	}
	ASSERT_EQ(tmpPosition, lastPosition);
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadRegionAtStartNonExactBeginingExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 1;
	graphite::position endPosition = 3728155;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 3728155;
	uint32_t totalCount = 12;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 1;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
		++count;
	}
	ASSERT_EQ(tmpPosition, lastPosition);
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadRegionAtStartExactBeginingExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 909434;
	graphite::position endPosition = 3728155;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 3728155;
	uint32_t totalCount = 12;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 1;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
		++count;
	}
	ASSERT_EQ(tmpPosition, lastPosition);
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadRegionAtEndExactBeginingExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 246020849;
	graphite::position endPosition = 248828840;
	graphite::position lastPosition = 248828840;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtEndNonExactBeginningExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 246020749;
	graphite::position endPosition = 248828840;
	graphite::position lastPosition = 248828840;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtEndExactBeginingNonExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 246020749;
	graphite::position endPosition = 300000000;
	graphite::position lastPosition = 248828840;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtEndNonExactBeginingNonExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition = 238828840;
	graphite::position endPosition = 348828840;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 248828840;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtMiddleExactBeginingExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition =  36314453;
	graphite::position endPosition = 38686990;
	graphite::position firstPosition =  36314453;
	graphite::position lastPosition = 38686990;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtMiddleNonExactBeginingExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition =  36314450;
	graphite::position endPosition = 38686990;
	graphite::position firstPosition =  36314453;
	graphite::position lastPosition = 38686990;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtMiddleExactBeginingNonExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition =  36314453;
	graphite::position endPosition = 38686991;
	graphite::position firstPosition =  36314453;
	graphite::position lastPosition = 38686990;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionAtMiddleNonExactBeginingNonExactEnd)
{
	std::string chrom = "1";
	graphite::position startPosition =  36314452;
	graphite::position endPosition = 38686991;
	graphite::position firstPosition =  36314453;
	graphite::position lastPosition = 38686990;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	graphite::position tmpPosition;
	variantListPtr->getNextVariant(variantPtr);
	ASSERT_EQ(variantPtr->getPosition(), firstPosition);
	while (variantListPtr->getNextVariant(variantPtr))
	{
		tmpPosition = variantPtr->getPosition(); // this will store the last position of the next variant, that way I have access to the last position
	}
	ASSERT_EQ(tmpPosition, lastPosition);
}

TEST(VCFFileReaderTests, ReadRegionNonExistantChrom)
{
	std::string chrom = "J";
	uint32_t totalCount = 0;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom;
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 0;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		++count;
	}
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadRegionWithNoVariants)
{
	std::string chrom = "1";
	graphite::position startPosition = 1;
	graphite::position endPosition = 100;
	graphite::position firstPosition = 909434;
	graphite::position lastPosition = 3728155;
	uint32_t totalCount = 0;
	std::string path = TEST_VCF_FILE;
    std::string regionString = chrom + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition);
	auto regionPtr = std::make_shared< graphite::Region >(regionString);
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	graphite::IVariant::SharedPtr variantPtr;
	uint32_t count = 0;
	graphite::position tmpPosition;
	while (variantListPtr->getNextVariant(variantPtr))
	{
		++count;
	}
	ASSERT_EQ(count, totalCount);
}

TEST(VCFFileReaderTests, ReadAllRegionsCheckCount)
{
	uint32_t totalCount = 8127;
	uint32_t count = 0;
	for (uint32_t i = 1; i <= 22; ++i)
	{
		std::string chrom = std::to_string(i);
		std::string path = TEST_VCF_FILE;
		std::string regionString = chrom;
		auto regionPtr = std::make_shared< graphite::Region >(regionString);
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(path, regionPtr);
		variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
		variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
		auto variantListPtr = variantManagerPtr->getCompleteVariantList();
		graphite::IVariant::SharedPtr variantPtr;
		graphite::position tmpPosition;
		while (variantListPtr->getNextVariant(variantPtr))
		{
			ASSERT_STREQ(chrom.c_str(), variantPtr->getChrom().c_str());
			++count;
		}
	}
	ASSERT_EQ(count, totalCount);
}

#endif // GRAPHITE_VCF_FILE_TESTS_HPP
