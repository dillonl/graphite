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

#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
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

		gwiz::vg::ReferenceNode::SharedPtr getReferenceVertexContainsPositionTest(gwiz::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

	};

	const std::string m_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(AdjudicatorTests, TestFindPositionSingleVertex)
	{
		gwiz::position pos = 6000;
		gwiz::position firstVariantPositionOffset = 10;
		gwiz::position secondVariantPositionOffset = 19;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		auto referenceNode = variantGraph->getReferenceVertexContainsPositionTest(6000);
		// ASSERT_EQ(true, true);
	}

	/*
	TEST(AdjudicatorTests, TestFindPosition)
	{
		gwiz::position pos = 6000;
		gwiz::position firstVariantPositionOffset = 10;
		gwiz::position secondVariantPositionOffset = 19;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		std::vector< std::string > variants;
		variants.push_back("G");
		variants.push_back("GA");
		testReferenceVariantGenerator.addVariant(pos + firstVariantPositionOffset, ".", 1, variants);
		testReferenceVariantGenerator.addVariant(pos + secondVariantPositionOffset, ".", 1, { "G" });

		auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
	}
	*/
}

#endif //GWIZ_ADJUDICATOR_VARIANT_GRAPH_TESTS_HPP
