#ifndef GWIZ_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
#define GWIZ_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
#include <chrono>
#include <thread>

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/adjudicator/graph/VariantContig.h"

#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"


namespace
{

	class VariantContigTest : public gwiz::adjudicator::VariantContig
	{
	public:
		typedef std::shared_ptr< VariantContigTest > VariantContigTestPtr;

		VariantContigTest(gwiz::Region::SharedPtr region, gwiz::vg::VariantGraph::GraphPtr graphPtr, gwiz::vg::VariantGraph::VariantVertexDescriptor startVertex, gwiz::vg::VariantGraph::VariantVertexDescriptor endVertex) :
			gwiz::adjudicator::VariantContig(region, graphPtr, startVertex, endVertex)
		{

		}

		~VariantContigTest()
		{
		}

		// gwiz::VariantGraph::Graph getGraph() { return this->m_graph_ptr; }

	};

	const std::string m_var_contig_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(VariantContigTests, TestBuildVariants)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 10;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_var_contig_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		std::vector< std::string > variants;
		variants.push_back("G");
		variants.push_back("GA");
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, variants);

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		auto referenceVertex = variantGraph->getReferenceVertexContainsPositionTest(6000);

		// auto region = std::make_shared< gwiz::Region >("20", pos, pos + m_var_contig_reference_string.size());
		// auto variantContig(region, variantGraph->getGraph(), )

		std::cout << "vertex desc: " << sizeof(referenceVertex) << std::endl;
		auto referenceNode = std::dynamic_pointer_cast< gwiz::vg::ReferenceNode >((*variantGraph->getGraph())[referenceVertex]);

		variantGraph->printGraph("test1.dot");

		/*
		bool refEq = memcmp(referenceNode->getSequence(), testReferenceVariantGenerator.getReference()->getSequence(), referenceNode->getLength())  == 0;
		ASSERT_TRUE(refEq);
		ASSERT_EQ(referenceNode->getPosition(), pos);
		*/
	}
}

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
