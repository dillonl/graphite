#ifndef GRAPHITE_ADJUDICATIONTESTS_HPP
#define GRAPHITE_ADJUDICATIONTESTS_HPP

#include "TestConfig.h"

#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariant.h"
#include "core/adjudicator/IAdjudicator.h"
#include "core/mapping/MappingManager.h"

#include "plugins/gssw/graph/GSSWGraph.h"
#include "plugins/gssw/graph/GSSWAdjudicator.h"
#include "plugins/gssw/graph/GSSWMapping.h"


#include <vector>

namespace
{
namespace adj_test
{
	using namespace graphite;
	using namespace graphite::gssw;

	class ReferenceTest : public IReference
	{
	public:
        ReferenceTest(Region::SharedPtr regionPtr, const char* sequence) : m_sequence(sequence) { this->m_region = regionPtr; }
		~ReferenceTest() {}

		const char* getSequence() override
		{
			return m_sequence;
		}

		size_t getSequenceSize() override
		{
			return strlen(m_sequence);
		}

	private:
		const char* m_sequence;
	};

	class VariantListTest : public IVariantList
	{
	public:
		VariantListTest(std::vector< IVariant::SharedPtr > variants) : m_variant_list(variants), m_index(0) {}
		~VariantListTest() {}

		bool getNextVariant(IVariant::SharedPtr& variantPtr) override
		{
			if (this->m_variant_list.size() <= this->m_index) { variantPtr = nullptr; return false; }
			else { variantPtr = this->m_variant_list[this->m_index++]; return true; }
		}

		size_t getCount() { return this->m_variant_list.size(); }
		void sort() {}
		void printToVCF(std::ostream& out) {}

	private:
		std::vector< IVariant::SharedPtr > m_variant_list;
		uint32_t m_index;
	};

	class AlignmentTest : public IAlignment
	{
	public:
		AlignmentTest(const std::string sequence, position pos) : m_sequence(sequence), m_position(pos) {}
		~AlignmentTest() {}

		const char* getSequence() override { return this->m_sequence.c_str(); }
		const size_t getLength() override { return this->m_sequence.size(); }
		const position getPosition() override { return this->m_position; }
	private:
		std::string m_sequence;
		position m_position;
	};

	GSSWGraph::SharedPtr getGSSWGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAdjudicator::SharedPtr adjudicatorPtr)
	{
		uint32_t startPosition = 1;
		uint32_t graphSize = 3000;
		auto gsswGraphPtr = std::make_shared< GSSWGraph >(referencePtr, variantListPtr, startPosition, graphSize, adjudicatorPtr->getMatchValue(), adjudicatorPtr->getMisMatchValue(), adjudicatorPtr->getGapOpenValue(), adjudicatorPtr->getGapExtensionValue());
		gsswGraphPtr->constructGraph();
		return gsswGraphPtr;
	}

	TEST(AdjudicationTest, AdjudicateDualSNPMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("A");
		auto altAllelePtr = std::make_shared< Allele >("G");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllelePtr };
		position pos = 10;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("CTCAAGTAGAATCTACTCTCTCAGGTGTTCATAATGTATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCTT", 2);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllelePtr->getTotalCount(), 1);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateDuoSNPMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("C");
		auto altAllelePtr = std::make_shared< Allele >("G");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllelePtr };
		position pos = 17;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		size_t alignmentOffset = 9;
		size_t alignmentSize = 150;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >(std::string(seq.c_str() + alignmentOffset, alignmentSize), alignmentOffset);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 1);
		ASSERT_EQ(altAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateTriSNPMisMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("T");
		auto altAllele1Ptr = std::make_shared< Allele >("G");
		auto altAllele2Ptr = std::make_shared< Allele >("A");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr, altAllele2Ptr };
		position pos = 15;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		size_t alignmentOffset = 9;
		size_t alignmentSize = 150;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >(std::string(seq.c_str() + alignmentOffset, alignmentSize), alignmentOffset);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 1);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 0);
		ASSERT_EQ(altAllele2Ptr->getTotalCount(), 0);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateTriSNPMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("T");
		auto altAllele1Ptr = std::make_shared< Allele >("G");
		auto altAllele2Ptr = std::make_shared< Allele >("A");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr, altAllele2Ptr };
		position pos = 15;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("AAATCAACTCTCTCAGGTGTTCATAATGTATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCTTATTCCTTC", 6);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 0);
		ASSERT_EQ(altAllele2Ptr->getTotalCount(), 1);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateShortDeletionMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");
		auto refAllelePtr = std::make_shared< Allele >("AGGTGTTCATAATGT");
		auto altAllele1Ptr = std::make_shared< Allele >("AT");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr };
		position pos = 24;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("TACTCTCTCATATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCTTATTCCTTCCCCTG", 6);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 1);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateShortInsertionMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");

		auto refAllelePtr = std::make_shared< Allele >("A");
		auto altAllele1Ptr = std::make_shared< Allele >("AGGGGGGGGGG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr };
		position pos = 5;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("ACTCAGGGGGGGGGGAGTAAAATCTACTCTCTCAGGTGTTCATAATGTATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCT", 1);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 1);
		ASSERT_EQ(gsswMappingPtr->getMappingScore(), (alignmentPtr->getLength() * match));
	}

	TEST(AdjudicationTest, AdjudicateShortInsertionMultipleMatch)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");

		auto refAllelePtr = std::make_shared< Allele >("A");
		auto altAllele1Ptr = std::make_shared< Allele >("AGGGGGGGGGG");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr };
		position pos = 5;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("ACTCAGGGGGGGGGGAGTAAAATCTACTCTCTCAGGTGTTCATAATGTATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCT", 1);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 1);
	}

	TEST(AdjudicationTest, AdjudicateDualIdenticalIntoVariant)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");

		auto refAllelePtr = std::make_shared< Allele >("TTCCCCTGACCCT");
		auto altAllele1Ptr = std::make_shared< Allele >("TTCCCCTGTTTTATCTTTT");
		auto altAllele2Ptr = std::make_shared< Allele >("TTCCGTTTTATCTTTT");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr, altAllele2Ptr };
		position pos = 157;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("GTTTTATCTTTTATTCTCCAAATTTCTTTCCAATTCATCTTTGTTCTTCCCTTTCCTTTTTACTCTCTTTAAACATTCTATGGACTCTGCCTCCTTCACACTGATATTGAACGCCCATAGTTTCATATTTTGGATTGCGATTGTTTTATT", 1);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 0);
		ASSERT_EQ(altAllele2Ptr->getTotalCount(), 0);
	}

	TEST(AdjudicationTest, AdjudicateDualIdenticalOutOfVariant)
	{
		std::string seq = TEST_REFERENCE_SEQUENCE;
		auto regionPtr = std::make_shared< Region >("1:1-3720");

		auto refAllelePtr = std::make_shared< Allele >("TGACCCTTCTTTTATTCTC");
		auto altAllele1Ptr = std::make_shared< Allele >("TGACCCTTACCTG");
		auto altAllele2Ptr = std::make_shared< Allele >("TGACCCTTGATTC");
		std::vector< IAllele::SharedPtr > altAllelePtrs = { altAllele1Ptr, altAllele2Ptr };
		position pos = 163;
		std::string chrom = "1";
		std::string dot = ".";
		auto variantPtr = std::make_shared< Variant >(pos, chrom, dot, dot, dot, refAllelePtr, altAllelePtrs);
		variantPtr->processOverlappingAlleles();
		std::vector< IVariant::SharedPtr > variantPtrs = { variantPtr };
		auto variantListPtr = std::make_shared< VariantList >(variantPtrs);
		auto referencePtr = std::make_shared< ReferenceTest >(regionPtr, seq.c_str());

		uint32_t percent = 80;
		int match = 1;
		int mismatch = 4;
		int gapOpen = 6;
		int gapExtension = 1;
		auto gsswAdjudicatorPtr = std::make_shared< GSSWAdjudicator >(percent, match, mismatch, gapOpen, gapExtension);
		auto gsswGraphPtr = getGSSWGraph(referencePtr, variantListPtr, gsswAdjudicatorPtr);
		auto alignmentPtr = std::make_shared< AlignmentTest >("CTCAGGTGTTCATAATGTATCAATGTATATTGCTTTAAGCCTGAAGGTAACCTAAGTAAAGATGTACCATGTTCCACCAATGCTTCTTTTGATCATCATTTTATCCTGTTTTTTCTTTAGGATTCTTTCTTATTCCTTCCCCTGACCCTT", 1);

		auto gsswMappingPtr = std::make_shared< GSSWMapping >(gsswGraphPtr->traceBackAlignment(alignmentPtr), alignmentPtr);
		gsswMappingPtr->printLongFormat();
		MappingManager::Instance()->registerMapping(gsswMappingPtr);
		MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicatorPtr);

		ASSERT_EQ(refAllelePtr->getTotalCount(), 0);
		ASSERT_EQ(altAllele1Ptr->getTotalCount(), 0);
		ASSERT_EQ(altAllele2Ptr->getTotalCount(), 0);
	}

}
}

#endif //GRAPHITE_ADJUDICATIONTESTS_HPP
