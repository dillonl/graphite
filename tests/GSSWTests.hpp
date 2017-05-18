#ifndef GRAPHITE_GSSWTESTS_HPP
#define GRAPHITE_GSSWTESTS_HPP

#include "TestClasses.hpp"

#include "gssw.h"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/reference/Reference.h"
#include "core/graph/GSSWGraph.h"
#include "core/file/FastaFileWriter.h"

#include <functional>

static uint32_t test_match = 1;
static uint32_t test_mismatch = 4;
static uint32_t gapOpen = 6;
static uint32_t gapExtension = 1;

void runTestsViaLambda(const std::string& alignment, const std::vector< gssw_node* >& nodes, int8_t* nt_table, int8_t* mat, std::function< void(gssw_graph_mapping* gm) > funct)
{
	gssw_graph* graph = gssw_graph_create(4);

	for (auto& node : nodes)
	{
		gssw_graph_add_node(graph, node);
	}

	gssw_graph_fill(graph, alignment.c_str(), alignment.size(), nt_table, mat, gapOpen, gapExtension, 15, 2);
	gssw_graph_mapping* gm = gssw_graph_trace_back (graph, alignment.c_str(), alignment.size(), test_match, test_mismatch, gapOpen, gapExtension);

	funct(gm);

	gssw_graph_mapping_destroy(gm);
	gssw_graph_destroy(graph);
	free(nt_table);
	free(mat);
}

TEST(GSSWTests, GSSWSimpleDiamondGraph)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(test_match, test_mismatch);

	std::string alignment = "ACGT";

	std::vector< std::string > nodeSeqs = {"A","G","C","GT"};
	std::vector< gssw_node* > nodes;
	uint32_t count = 0;
	for (auto& nodeSeq : nodeSeqs)
	{
		nodes.emplace_back(gssw_node_create(nullptr, count++, nodeSeq.c_str(), nt_table, mat));
	}

	gssw_nodes_add_edge(nodes[0], nodes[1]);
	gssw_nodes_add_edge(nodes[0], nodes[2]);
	gssw_nodes_add_edge(nodes[1], nodes[3]);
	gssw_nodes_add_edge(nodes[2], nodes[3]);

	auto funct = [](gssw_graph_mapping* gm)
	{
		std::vector< uint32_t > path = {0,2,3};
		gssw_node_cigar* nc = gm->cigar.elements;
		for (int i = 0; i < gm->cigar.length; ++i, ++nc)
		{
			ASSERT_EQ(path[i], nc->node->id);
		}
	};
	runTestsViaLambda(alignment, nodes, nt_table, mat, funct);
}

TEST(GSSWTests, GSSWSimple3WayGraph)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(test_match, test_mismatch);

	std::string alignment = "ACGT";

	std::vector< std::string > nodeSeqs = {"A","G","T","C","GT"};
	std::vector< gssw_node* > nodes;
	uint32_t count = 0;
	for (auto& nodeSeq : nodeSeqs)
	{
		nodes.emplace_back(gssw_node_create(nullptr, count++, nodeSeq.c_str(), nt_table, mat));
	}

	gssw_nodes_add_edge(nodes[0], nodes[1]);
	gssw_nodes_add_edge(nodes[0], nodes[2]);
	gssw_nodes_add_edge(nodes[0], nodes[3]);
	gssw_nodes_add_edge(nodes[1], nodes[4]);
	gssw_nodes_add_edge(nodes[2], nodes[4]);
	gssw_nodes_add_edge(nodes[3], nodes[4]);

	auto funct = [](gssw_graph_mapping* gm)
		{
			std::vector< uint32_t > path = {0,3,4};
			gssw_node_cigar* nc = gm->cigar.elements;
			for (int i = 0; i < gm->cigar.length; ++i, ++nc)
			{
				ASSERT_EQ(path[i], nc->node->id);
			}
		};
	runTestsViaLambda(alignment, nodes, nt_table, mat, funct);
}

TEST(GSSWTests, GSSWComplexDiamondGraph)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(test_match, test_mismatch);

	std::string alignment = "GGATT";

	std::vector< std::string > nodeSeqs = {"G","A","G","GATT"};
	std::vector< gssw_node* > nodes;
	uint32_t count = 0;
	for (auto& nodeSeq : nodeSeqs)
	{
		nodes.emplace_back(gssw_node_create(nullptr, count++, nodeSeq.c_str(), nt_table, mat));
	}

	gssw_nodes_add_edge(nodes[0], nodes[1]);
	gssw_nodes_add_edge(nodes[0], nodes[2]);
	gssw_nodes_add_edge(nodes[1], nodes[3]);
	gssw_nodes_add_edge(nodes[2], nodes[3]);

	auto funct = [](gssw_graph_mapping* gm)
		{
			std::vector< uint32_t > path = {2,3};
			gssw_node_cigar* nc = gm->cigar.elements;
			for (int i = 0; i < gm->cigar.length; ++i, ++nc)
			{
				ASSERT_EQ(path[i], nc->node->id);
			}
		};
	runTestsViaLambda(alignment, nodes, nt_table, mat, funct);
}

TEST(GSSWTests, GSSWSingleton)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(test_match, test_mismatch);

	std::string alignment = "AAAAAAAGAAAA";

	std::vector< std::string > nodeSeqs = {"TAAAAAAAGAAAAT"};
	std::vector< gssw_node* > nodes;
	uint32_t count = 0;
	for (auto& nodeSeq : nodeSeqs)
	{
		nodes.emplace_back(gssw_node_create(nullptr, count++, nodeSeq.c_str(), nt_table, mat));
	}

	auto funct = [](gssw_graph_mapping* gm)
		{
			gssw_node_cigar* nc = gm->cigar.elements;
			for (int i = 0; i < gm->cigar.length; ++i, ++nc)
			{
				for (int j = 0; j < nc->cigar->length; ++j)
				{
					ASSERT_EQ(nc->cigar->elements[i].length, 12);
					ASSERT_EQ(nc->cigar->elements[i].type, 'M');
				}
			}
		};
	runTestsViaLambda(alignment, nodes, nt_table, mat, funct);
}

/*
TEST(GSSWTests, FastaFileWriter)
{
    std::vector< std::string > headers {">chr1:1000:0", ">chr1:1000:1"};
    std::vector< std::string > sequences {"ATCAT", "GGGCC"};

    std::cout << "FastaFileWriter test is running!" << std::endl;
    graphite::FastaFileWriter fastaFileWriter;

    std::cout << "File opened: " << fastaFileWriter.open("gtestFastaFile.fa") << "\n";
    fastaFileWriter.write(headers, sequences);
    fastaFileWriter.close();
    std::cout << "File closed: " << "\n";
}
*/

// Need to finish adding the TestReferenc1 class.
// Briefly look at the IReference interface.
// Finish updating the SimplePaths test.

class TestReference1: public graphite::IReference
{
private:

public:
    TestReference1 () {}

    ~TestReference1 () {}

    const char* getSequence() override { return m_sequence.c_str(); }
    size_t getSequenceSize() override { return m_sequence.size(); }

    void setSequence (graphite::Region::SharedPtr regionPtr, const std::string& sequence)
    {
        m_region = regionPtr;
        m_sequence = sequence;
    }
};

// The problem lies in GSSWGraph->constructGraph() with the referenceSequenceString line.
//   In the getSequenceFromRegion function underflow occurs in the if statement when 1 is subtracted from startPosition.
//   The startPosition for both regionPtr and referencePtr = 1. This causes underflow when execution enters the if statement checking for BASED.
// Do a traceback starting with this line.
// Notes:
//   currentReferencePosition is set equal to regionPtr startPosition before entering whiel loop (GSSWGraph.cpp).
// Should currentReferencePosition = 0 ?
//
// The solution likely lies with adjusting the parameters of either regionPtr or referencePtr.
//   Why doesn't currentReferencePosition change when I change the startPosition of regionPtr from 1 to 2?
//   Try printing out all local variables in (gdb) info locals.
//   Can compile in gdb using make.

TEST(GSSWTests, GSSWSimplePaths)
{
    std::string VCF_LINE_1 = "chr1\t10\trs11575897\tA\tT\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line.
    auto regionPtr = std::make_shared< graphite::Region > ("1", 1, 50, graphite::Region::BASED::ONE);
    std::string sequence = "TGAAGGCCAAAATTCAGATTCAGGACCCCTCCCGGGTAAAAATATATATA";

    auto referencePtr = std::make_shared< TestReference1 > ();
    referencePtr->setSequence(regionPtr, sequence);

    /*
    std::cout << "Start position from region: " << regionPtr->getStartPosition() << std::endl;
    std::cout << "Start position from reference: " << referencePtr->getRegion()->getStartPosition() << std::endl;
    std::cout << "getRegionString: " << regionPtr->getRegionString() << std::endl;
    std::cout << "getEndPosition(): " << regionPtr->getEndPosition() << std::endl;
    std::cout << "getBased(): " << static_cast<int>(regionPtr->getBased()) << std::endl;
    std::cout << "getBased() from referencePtr: " << static_cast<int>(referencePtr->getRegion()->getBased()) << std::endl;
    std::cout << "region string from referencePtr: " << referencePtr->getRegion()->getRegionString() << std::endl;
    //std::cout << "getSequence test: " << std::endl;
    // The following test doesn't work.
    //std::string referenceString = referencePtr->getSequenceFromRegion(regionPtr);
    //std::cout << "The reference string: \n" << sequence << std::endl;
    */
    uint32_t readLength = 5;  // Need to double check that this works with BuildVariant
    auto variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), referencePtr, readLength);

    std::vector< graphite::IVariant::SharedPtr > variantPtrs = {variantPtr};
    auto variantListPtr = std::make_shared< graphite::VariantList >(variantPtrs, referencePtr);
    uint32_t numGraphCopies = 1;    // Dummy variable for the updated GSSWGraph constructor.
    auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, regionPtr, 1, 1, 1, 1, numGraphCopies);

    std::cout << "Created gsswGraphPtr" << std::endl;
    gsswGraphPtr->constructGraph();
    
    /*
     * Segmentation fault occurs when the constructGraph function is used.
     * constructGraph uses region, reference and variant ptrs to get information.
     * I think regionPtr is setup correctly.
     * I think referencePtr is setup correctly.
     * I think variantPtr is setup correctly.
     * I think variantListPtr is setup correctly.
    std::cout << "Constructed graph" << std::endl;
    std::vector< std::string > graphPathHeaders = gsswGraphPtr->getGraphPathHeaders();
    std::vector< std::string > graphPathSequences = gsswGraphPtr->getGraphPathSequences();

    std::cout << "Number of headers: " << graphPathHeaders.size() << std::endl;
    std::cout << "Number of sequences: " << graphPathSequences.size() << std::endl;
    for (int i = 0; i < graphPathHeaders.size() && i < graphPathSequences.size(); ++i)
    {
        std::cout << graphPathHeaders[i] << std::endl;
        std::cout << graphPathSequences[i] << std::endl;
    }
    */
}

/*
TEST(GSSWTests, TestPaths)
{
	std::vector< std::string > v1 = {"1"};
	std::vector< std::string > v2 = {"2a", "2b", "2c", "2d"};
	std::vector< std::string > v3 = {"3"};
	std::vector< std::string > v4 = {"4a", "4b", "4c", "4d"};
	std::vector< std::string > v5 = {"5"};
	std::vector< std::vector< std::string > > pathsComponents = {v1, v2, v3, v4, v5};
	std::vector< std::string > paths;
	graphite::GSSWGraph::generatePaths(pathsComponents, paths, 0, "");
	for (auto path : paths)
	{
		std::cout << path << std::endl;
	}

}
*/

#endif //GRAPHITE_GSSWTESTS_HPP
