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

	TEST(VariantContigTests, TestBuildVariantsMultiSingleSite)
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

		gwiz::position startPos = 6001;
		gwiz::position endPos = 6025;

		auto startVertex = variantGraph->getReferenceVertexContainsPositionTest(startPos);
		auto endVertex = variantGraph->getReferenceVertexContainsPositionTest(endPos);

		std::string regionString = "20:" + to_string(startPos) + "-" + to_string(endPos);
		auto region = std::make_shared< gwiz::Region >(regionString);
		auto contig = std::make_shared< VariantContigTest >(region, variantGraph->getGraph(), startVertex, endVertex);
		contig->buildVariantContig();
		auto contigs = contig->getContigs();
		std::vector< std::string > contigStrings;
		for_each (contigs.begin(), contigs.end(), [&](std::tuple< std::list< gwiz::vg::VariantGraph::VariantVertexDescriptor >, std::string >& con) {
				contigStrings.push_back(std::get< 1 >(con));
			});
		ASSERT_EQ(contigs.size(), 3);
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTTATCCACCCACCTTG") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTGAATCCACCCACCTTG") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTGATCCACCCACCTTG") != contigStrings.end());
	}

	TEST(VariantContigTests, TestBuildVariantsMultipleSite)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 10;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_var_contig_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1,{"A"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset + 1, "2", 1, {"G"});

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();

		gwiz::position startPos = 6001;
		gwiz::position endPos = 6025;

		auto startVertex = variantGraph->getReferenceVertexContainsPositionTest(startPos);
		auto endVertex = variantGraph->getReferenceVertexContainsPositionTest(endPos);

		std::string regionString = "20:" + to_string(startPos) + "-" + to_string(endPos);
		auto region = std::make_shared< gwiz::Region >(regionString);
		auto contig = std::make_shared< VariantContigTest >(region, variantGraph->getGraph(), startVertex, endVertex);
		contig->buildVariantContig();
		auto contigs = contig->getContigs();
		std::vector< std::string > contigStrings;
		for_each (contigs.begin(), contigs.end(), [&](std::tuple< std::list< gwiz::vg::VariantGraph::VariantVertexDescriptor >, std::string >& con) {
				contigStrings.push_back(std::get< 1 >(con));
			});
		ASSERT_EQ(contigs.size(), 4);
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTTATCCACCCACCTTG") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTTGTCCACCCACCTTG") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTAATCCACCCACCTTG") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AATAAAAGTTAGTCCACCCACCTTG") != contigStrings.end());
	}
}

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
