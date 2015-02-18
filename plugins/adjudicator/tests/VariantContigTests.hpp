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

	};

	const std::string m_var_contig_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(VariantContigTests, TestBuildVariants)
	{
		gwiz::position pos = 6000;
		gwiz::position firstVariantPositionOffset = 10;
		gwiz::position secondVariantPositionOffset = 19;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_var_contig_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		auto referenceVertex = variantGraph->getReferenceVertexContainsPositionTest(6000);
		std::cout << "vertex desc: " << sizeof(referenceVertex) << std::endl;
		auto referenceNode = std::dynamic_pointer_cast< gwiz::vg::ReferenceNode >((*variantGraph->getGraph())[referenceVertex]);

		bool refEq = memcmp(referenceNode->getSequence(), testReferenceVariantGenerator.getReference()->getSequence(), referenceNode->getLength())  == 0;
		ASSERT_TRUE(refEq);
		ASSERT_EQ(referenceNode->getPosition(), pos);
	}
}

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
