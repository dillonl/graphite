#ifndef GRAPHITE_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
#define GRAPHITE_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
#include <chrono>
#include <thread>

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/adjudicator/graph/VariantContig.h"

#include "core/variant/VCFFileReader.h"
#include "core/variant/IVariant.h"
#include "core/reference/FastaReference.h"


namespace
{

	class VariantContigTest : public graphite::adjudicator::VariantContig
	{
	public:
		typedef std::shared_ptr< VariantContigTest > VariantContigTestPtr;

		VariantContigTest(uint32_t padding, graphite::vg::VariantGraph::GraphPtr graphPtr, graphite::vg::VariantGraph::VariantVertexDescriptor startVertex, graphite::vg::VariantGraph::VariantVertexDescriptor endVertex) :
			graphite::adjudicator::VariantContig(padding, graphPtr, startVertex, endVertex)
		{

		}

		~VariantContigTest()
		{
		}

		// graphite::VariantGraph::Graph getGraph() { return this->m_graph_ptr; }

	};

	const std::string m_var_contig_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(VariantContigTests, TestBuildVariantsMultiSingleSite)
	{
		graphite::position pos = 6000;
		graphite::position variantPositionOffset = 10;
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_var_contig_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		std::vector< std::string > variants;
		variants.push_back("G");
		variants.push_back("GA");
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, variants);

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();

		graphite::position startPos = 6001;
		graphite::position endPos = 6025;

		auto startVertex = variantGraph->getReferenceVertexContainsPositionTest(startPos);
		auto endVertex = variantGraph->getReferenceVertexContainsPositionTest(endPos);

		auto contig = std::make_shared< VariantContigTest >(5, variantGraph->getGraph(), startVertex, endVertex);
		contig->buildVariantContig();
		auto contigs = contig->getContigs();
		std::vector< std::string > contigStrings;
		for_each (contigs.begin(), contigs.end(), [&](const graphite::adjudicator::VariantContig::ContigTuplePtr& con) {
				contigStrings.push_back(std::get< 1 >(*con));
			});
		ASSERT_EQ(contigs.size(), 3);
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTTATCCA") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTGAATCCA") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTGATCCA") != contigStrings.end());

	}

	TEST(VariantContigTests, TestBuildVariantsMultipleSite)
	{
		graphite::position pos = 6000;
		graphite::position variantPositionOffset = 10;
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_var_contig_reference_string, "20", pos);
		//AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1,{"A"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset + 1, "2", 1, {"G"});

		auto variantGraph = std::make_shared< AdjudicatorGraphTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();

		graphite::position startPos = 6001;
		graphite::position endPos = 6025;

		auto startVertex = variantGraph->getReferenceVertexContainsPositionTest(startPos);
		auto endVertex = variantGraph->getReferenceVertexContainsPositionTest(endPos);

		// std::string regionString = "20:" + to_string(startPos) + "-" + to_string(endPos);
		// auto region = std::make_shared< graphite::Region >(regionString);

		auto contig = std::make_shared< VariantContigTest >(5, variantGraph->getGraph(), startVertex, endVertex);
		contig->buildVariantContig();
		auto contigs = contig->getContigs();
		std::vector< std::string > contigStrings;
		for_each (contigs.begin(), contigs.end(), [&](const graphite::adjudicator::VariantContig::ContigTuplePtr& con) {
				contigStrings.push_back(std::get< 1 >(*con));
			});
		ASSERT_EQ(contigs.size(), 4);
		//      AAAGTTATCCAC
		//      AAAGTTGTCCAC
		//		AAAGTAATCCAC
		//		AAAGTAGTCCAC
		// AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTTATCCAC") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTTGTCCAC") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTAATCCAC") != contigStrings.end());
		ASSERT_TRUE(std::find(contigStrings.begin(), contigStrings.end(), "AAAGTAGTCCAC") != contigStrings.end());
	}
}

#endif //GRAPHITE_ADJUDICATOR_VARIANTCONTIG_TESTS_HPP
