#ifndef GRAPHITE_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
#define GRAPHITE_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
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

	class AdjudicatorGraphTest : public graphite::adjudicator::AdjudicatorGraph
	{
	public:
		typedef std::shared_ptr< AdjudicatorGraphTest > AdjudicatorGraphTestPtr;

		AdjudicatorGraphTest(graphite::IReference::SharedPtr referencePtr, graphite::IVariantList::SharedPtr variantListPtr) :
			graphite::adjudicator::AdjudicatorGraph(referencePtr, variantListPtr)
		{
		}

		~AdjudicatorGraphTest()
		{
		}

		graphite::vg::VariantGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(graphite::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		graphite::vg::VariantGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}

	};

	const std::string m_adj_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	/*
	TEST(AdjudicatorTests, TestFindPositionSingleVertex)
	{
		graphite::position pos = 6000;
		graphite::position firstVariantPositionOffset = 10;
		graphite::position secondVariantPositionOffset = 19;
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_adj_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		auto referenceVertex = variantGraph->getReferenceVertexContainsPositionTest(6000);
		auto referenceNode = std::dynamic_pointer_cast< graphite::vg::ReferenceNode >((*variantGraph->getGraph())[referenceVertex]);

		bool refEq = memcmp(referenceNode->getSequence(), testReferenceVariantGenerator.getReference()->getSequence(), referenceNode->getLength())  == 0;
		ASSERT_TRUE(refEq);
		ASSERT_EQ(referenceNode->getPosition(), pos);
	}

	TEST(AdjudicatorTests, TestFindPosition)
	{
		graphite::position pos = 6000;
		graphite::position firstVariantPositionOffset = 10;
		graphite::position secondVariantPositionOffset = 19;
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_adj_reference_string, "20", pos);
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


		graphite::Region::SharedPtr regionPtr = std::make_shared< graphite::Region >("20");

		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;

		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReader = std::make_shared<graphite::VCFFileReader>(vcfPath);
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


		// graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< graphite::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}
}

#endif //GRAPHITE_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
