#ifndef GWIZ_GSSWTESTS_HPP
#define GWIZ_GSSWTESTS_HPP

#include "gssw/gssw.h"

#include <functional>

static uint32_t match = 1;
static uint32_t mismatch = 4;
static uint32_t gapOpen = 6;
static uint32_t gapExtension = 1;

void runTestsViaLambda(const std::string& alignment, const std::vector< gssw_node* >& nodes, int8_t* nt_table, int8_t* mat, std::function< void(gssw_graph_mapping* gm) > funct)
{
	gssw_graph* graph = gssw_graph_create(4);

	for (auto& node : nodes)
	{
		gssw_graph_add_node(graph, node);
	}

	gssw_graph_fill(graph, alignment.c_str(), nt_table, mat, gapOpen, gapExtension, 15, 2);
	gssw_graph_mapping* gm = gssw_graph_trace_back (graph, alignment.c_str(), alignment.size(), match, mismatch, gapOpen, gapExtension);

	funct(gm);

	gssw_graph_mapping_destroy(gm);
	gssw_graph_destroy(graph);
	free(nt_table);
	free(mat);
}

TEST(GSSWTests, GSSWSimpleDiamondGraph)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(match, mismatch);

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
	int8_t* mat = gssw_create_score_matrix(match, mismatch);

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
	int8_t* mat = gssw_create_score_matrix(match, mismatch);

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

/*
//current gssw fails this test

TEST(GSSWTests, GSSWComplexDiamondGraph2)
{
	int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(match, mismatch);

	std::string alignment = "ATGCTCAATGCGTAG";
	std::vector< std::string > nodeSeqs = {"ATGCTCG","A","C","ATGCGTAG"};
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
			std::vector< uint32_t > path = {0,1,3};
			std::string cigar = "";
			std::string cigarMapping = "7M1M8M";
			gssw_node_cigar* nc = gm->cigar.elements;
			for (int i = 0; i < gm->cigar.length; ++i, ++nc)
			{
				for (int j = 0; j < nc->cigar->length; ++j)
				{
					cigar += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
				}
				ASSERT_EQ(path[i], nc->node->id);
			}
			ASSERT_STREQ(cigarMapping.c_str(), cigar.c_str());
		};
	runTestsViaLambda(alignment, nodes, nt_table, mat, funct);
}
*/

#endif //GWIZ_GSSWTESTS_HPP
