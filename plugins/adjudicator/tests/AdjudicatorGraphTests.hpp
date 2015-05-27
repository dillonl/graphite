#ifndef GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
#define GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
#include <chrono>
#include <thread>

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/adjudicator/graph/AdjudicatorGraph.h"

#include "core/variant/VCFFileReader.h"
#include "core/variant/IVariant.h"
#include "core/reference/FastaReference.h"


namespace
{

	class AdjudicatorGraphTest : public gwiz::adjudicator::AdjudicatorGraph
	{
	public:
		typedef std::shared_ptr< AdjudicatorGraphTest > AdjudicatorGraphTestPtr;

		AdjudicatorGraphTest(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr) :
			gwiz::adjudicator::AdjudicatorGraph(referencePtr, variantListPtr)
		{
		}

		~AdjudicatorGraphTest()
		{
		}

		gwiz::vg::VariantGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(gwiz::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		gwiz::vg::VariantGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}

	};

	const std::string m_adj_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	/*
	TEST(AdjudicatorTests, TestFindPositionSingleVertex)
	{
		gwiz::position pos = 6000;
		gwiz::position firstVariantPositionOffset = 10;
		gwiz::position secondVariantPositionOffset = 19;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_adj_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		auto referenceVertex = variantGraph->getReferenceVertexContainsPositionTest(6000);
		auto referenceNode = std::dynamic_pointer_cast< gwiz::vg::ReferenceNode >((*variantGraph->getGraph())[referenceVertex]);

		bool refEq = memcmp(referenceNode->getSequence(), testReferenceVariantGenerator.getReference()->getSequence(), referenceNode->getLength())  == 0;
		ASSERT_TRUE(refEq);
		ASSERT_EQ(referenceNode->getPosition(), pos);
	}

	TEST(AdjudicatorTests, TestFindPosition)
	{
		gwiz::position pos = 6000;
		gwiz::position firstVariantPositionOffset = 10;
		gwiz::position secondVariantPositionOffset = 19;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_adj_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		std::vector< std::string > variants;
		variants.push_back("G");
		variants.push_back("GA");
		testReferenceVariantGenerator.addVariant(pos + firstVariantPositionOffset, ".", 1, variants);
		testReferenceVariantGenerator.addVariant(pos + secondVariantPositionOffset, ".", 1, { "G" });

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
	}

	*/

	TEST(AdjudicatorTests, TestConstructChr20)
	{


		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20");

		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;

		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(vcfPath);
		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(fastaReferencePtr, vcfFileReader);
		variantGraph->constructGraph();

		std::cout << "contigs constructed" << std::endl;

		auto contigs = variantGraph->getContigs();
		size_t contigCount = 0;
		while (!contigs.empty())
		{
			auto contig = contigs.front();
			contigCount += contig->getContigs().size();
			contigs.pop();
		}
		std::cout << "Contig Count: " << contigCount << std::endl;

		// variantGraph->printGraph("test1.dot");


		// gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}
}

#endif //GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
