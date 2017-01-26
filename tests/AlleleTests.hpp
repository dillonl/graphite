#ifndef GRAPHITE_ALLELETESTS_HPP
#define GRAPHITE_ALLELETESTS_HPP

#include "TestConfig.h"

#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariant.h"
#include "core/adjudicator/IAdjudicator.h"
#include "core/mapping/MappingManager.h"

#include "core/graph/GSSWGraph.h"
#include "core/adjudicator/GSSWAdjudicator.h"
#include "core/mapping/GSSWMapping.h"


#include <vector>

namespace
{
namespace adj_test
{
	using namespace graphite;

	TEST(AlleleTest, overlappingPrefixCountZero)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("TGACCCTTCTTTTATTCTC");
		auto altAllelePtr = std::make_shared< Allele >("ACCTG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllelePtr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllelePtr), 0);
	}

	TEST(AlleleTest, overlappingPrefixCountNoSuffix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("TGACCCTTCTTTTATTCTCT");
		auto altAllele1Ptr = std::make_shared< Allele >("TGACCCTTACCTG");
		auto altAllele2Ptr = std::make_shared< Allele >("TGACCCTTGATTCA");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 8);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 8);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 8);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 0);
	}

	TEST(AlleleTest, overlappingSuffixCountNoPrefix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("ATCGACGG");
		auto altAllele1Ptr = std::make_shared< Allele >("GTCG");
		auto altAllele2Ptr = std::make_shared< Allele >("TAAAAAAAG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 1);
	}

	TEST(AlleleTest, overlappingSuffixCountPrefix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("GATCGACGG");
		auto altAllele1Ptr = std::make_shared< Allele >("GTCG");
		auto altAllele2Ptr = std::make_shared< Allele >("GCAAAAAAAG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 1);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 1);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 1);
	}

	TEST(AlleleTest, overlappingMultiMatchOneMisMatchCountPrefix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("ATCGACGG");
		auto altAllele1Ptr = std::make_shared< Allele >("GTC");
		auto altAllele2Ptr = std::make_shared< Allele >("GCAAAAAAA");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 1);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 0);
	}

	TEST(AlleleTest, overlappingMultiMatchOneMisMatchCountSuffix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("ATCGACGGA");
		auto altAllele1Ptr = std::make_shared< Allele >("TC");
		auto altAllele2Ptr = std::make_shared< Allele >("GCAAAAAAA");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 1);
	}

	TEST(AlleleTest, overlappingMultiMatchPrefixMixedCountSuffix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("TCGACGGA");
		auto altAllele1Ptr = std::make_shared< Allele >("TC");
		auto altAllele2Ptr = std::make_shared< Allele >("TAAAAAAA");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 2);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 2);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 1);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 1);
	}

	TEST(AlleleTest, overlappingMultiMatchPrefixMixedHigherCountOnPrefix)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("TCGACGGA");
		auto altAllele1Ptr = std::make_shared< Allele >("TCTTTTT");
		auto altAllele2Ptr = std::make_shared< Allele >("TCGACAAAAAG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr,altAllele2Ptr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();

		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(refAllelePtr), 5);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele1Ptr), 2);
		ASSERT_EQ(variantPtr->getAllelePrefixOverlapMaxCount(altAllele2Ptr), 5);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(refAllelePtr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele1Ptr), 0);
		ASSERT_EQ(variantPtr->getAlleleSuffixOverlapMaxCount(altAllele2Ptr), 0);
	}
}
}
#endif //GRAPHITE_ALLELETESTS_HPP
