#ifndef GRAPHITE_FASTAREFERENCETESTS_HPP
#define GRAPHITE_FASTAREFERENCETESTS_HPP

#include "TestConfig.h"

#include "core/region/Region.h"
#include "core/reference/FastaReference.h"

TEST(FastaTests, GetBaseOneRegion)
{
	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);
	auto regionPtr = std::make_shared< graphite::Region >("1:1-10", graphite::Region::BASED::ONE);
	std::string regionString = referencePtr->getSequenceFromRegion(regionPtr);
	ASSERT_STREQ(regionString.c_str(), "CTATGATGTT");
}

TEST(FastaTests, GetBaseOneRegionSimple)
{
	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);
	auto regionPtr = std::make_shared< graphite::Region >("1:2-3", graphite::Region::BASED::ONE);
	std::string regionString = referencePtr->getSequenceFromRegion(regionPtr);
	ASSERT_STREQ(regionString.c_str(), "TA");
}

TEST(FastaTests, GetBaseZeroRegion)
{
	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ZERO);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);
	auto regionPtr = std::make_shared< graphite::Region >("1:0-10", graphite::Region::BASED::ZERO);
	std::string regionString = referencePtr->getSequenceFromRegion(regionPtr);
	ASSERT_STREQ(regionString.c_str(), "CTATGATGTT");
}

TEST(FastaTests, GetBaseZeroRegionSimple)
{
	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ZERO);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);
	auto regionPtr = std::make_shared< graphite::Region >("1:2-3", graphite::Region::BASED::ZERO);
	std::string regionString = referencePtr->getSequenceFromRegion(regionPtr);
	ASSERT_STREQ(regionString.c_str(), "A");
}

#endif
