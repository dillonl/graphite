#ifndef GRAPHITE_GSSWTESTS_HPP
#define GRAPHITE_GSSWTESTS_HPP

#include "TestClasses.hpp"

#include "gssw.h"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/reference/Reference.h"
#include "core/graph/GSSWGraph.h"

#include <functional>

static uint32_t test_match = 1;
static uint32_t test_mismatch = 4;
static uint32_t gapOpen = 6;
static uint32_t gapExtension = 1;

/*
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
*/

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
