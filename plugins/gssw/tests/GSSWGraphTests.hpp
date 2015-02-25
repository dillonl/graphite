#ifndef GWIZ_GSSW_GSSWGRAPH_TESTS_HPP
#define GWIZ_GSSW_GSSWGRAPH_TESTS_HPP
#include <chrono>
#include <thread>

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/gssw/graph/GSSWGraph.h"

#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"


namespace
{

	class GSSWGraphTest : public gwiz::gssw::GSSWGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraphTest > GSSWGraphTestPtr;

		GSSWGraphTest(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr) :
			gwiz::gssw::GSSWGraph(referencePtr, variantListPtr)
		{
		}

		~GSSWGraphTest()
		{
		}

		/*
		gwiz::gssw::GSSWGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(gwiz::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		gwiz::gssw::GSSWGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}
		*/

	};

	const std::string m_gssw_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(GSSWTests, TestConstructChr20)
	{
		ASSERT_TRUE(true);

		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20");

		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;

		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(vcfPath);
		auto gsswGraph = std::make_shared< GSSWGraphTest >(fastaReferencePtr, vcfFileReader);
		gsswGraph->constructGraph();
/*
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
*/

		// variantGraph->printGraph("test1.dot");


		// gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}

}

#endif //GWIZ_GSSW_GSSWGRAPH_TESTS_HPP
