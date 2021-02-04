#pragma once

#include "Graph.h"
#include "Node.h"
#include "core/util/Noncopyable.hpp"

#include "gssw.h"

namespace graphite
{
	class GraphTraceback : Noncopyable
	{
	public:
		GraphTraceback(Graph::SharedPtr graphPtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue) : m_graph_ptr(graphPtr), m_match_value(matchValue), m_mismatch_value(mismatchValue), m_gap_open_value(gapOpenValue), m_gap_extension_value(gapExtensionValue)
		{
		}

		~GraphTraceback()
		{
		}

		void processGraph(Alignment::SharedPtr alignmentPtr)
		{
			auto nodePtrsMap = m_graph_ptr->getNodePtrsMap();
			gssw_sse2_disable();
			int8_t* nt_table = gssw_create_nt_table();
			int8_t* mat = gssw_create_score_matrix(m_match_value, m_mismatch_value);
			gssw_graph* graph = gssw_graph_create(nodePtrsMap.size());

			unordered_map< uint32_t, gssw_node* > gsswNodePtrsMap;
			Node::SharedPtr nodePtr = m_graph_ptr->getFirstNode();

			// std::lock_guard< std::mutex > l(m_graph_mutex);
			// std::cout << "node from graph: " << nodePtr->getAllelePtr() << std::endl;
			gssw_node* gsswNode = (gssw_node*)gssw_node_create(nodePtr.get(), nodePtr->getID(), nodePtr->getSequence().c_str(), nt_table, mat);
			gsswNodePtrsMap.emplace(nodePtr->getID(), gsswNode);
			gssw_graph_add_node(graph, gsswNode);
			while (nodePtr != nullptr)
			{
				Node::SharedPtr nextRefNodePtr = nullptr;
				for (auto outNodePtr : nodePtr->getOutNodes())
				{
					if (nodePtrsMap.find(outNodePtr->getID()) == nodePtrsMap.end())
					{
						continue; // skip removed nodes
					}
					if (outNodePtr->getAlleleType() == Node::ALLELE_TYPE::REF)
					{
						nextRefNodePtr = outNodePtr;
					}
					gsswNode = (gssw_node*)gssw_node_create(outNodePtr.get(), outNodePtr->getID(), outNodePtr->getSequence().c_str(), nt_table, mat);
					gssw_graph_add_node(graph, gsswNode);
					gsswNodePtrsMap.emplace(outNodePtr->getID(), gsswNode);
				}
				nodePtr = nextRefNodePtr;
			}

			for (auto iter : nodePtrsMap)
			{
				Node::SharedPtr nodePtr = iter.second;
				auto gsswIter = gsswNodePtrsMap.find(nodePtr->getID());
				if (gsswIter == gsswNodePtrsMap.end())
				{
					continue;
				}
				gssw_node* gsswNode = gsswIter->second;
				for (auto outNodePtr : nodePtr->getOutNodes())
				{
					auto iter = gsswNodePtrsMap.find(outNodePtr->getID());
					if (iter != gsswNodePtrsMap.end())
					{
						gssw_node* gsswOutNode = iter->second;
						gssw_nodes_add_edge(gsswNode, gsswOutNode);
					}
				}
			}
			gssw_graph_fill(graph, alignmentPtr->getSequence(), nt_table, mat, m_gap_open_value, m_gap_extension_value, 0, 0, 15, 2, true);
			gssw_graph_mapping* gm = gssw_graph_trace_back (graph, alignmentPtr->getSequence(), alignmentPtr->getLength(), nt_table, mat, m_gap_open_value, m_gap_extension_value, 0, 0);
			processTraceback(gm, alignmentPtr->getLength());

			setNormalizedCigarString(gm);

			gssw_graph_mapping_destroy(gm);

			// note that nodes which are referred to in this graph are destroyed as well
			gssw_graph_destroy(graph);

			free(nt_table);
			free(mat);
		}

		std::vector< Node* > getTracebackNodePtrs() { return m_traceback_node_ptrs; }
		uint32_t getTotalScore() { return m_total_score; }
		uint32_t getSoftClipOccurrences() { return m_soft_clip_occurrences; }
		std::string getCigarString() { return m_cigar_string; }
		uint32_t getNodeScorePercent(Node* nodePtr)
		{
			auto iter = m_node_id_node_score_map.find(nodePtr->getID());
			if (iter == m_node_id_node_score_map.end()) { return -1; }
			return iter->second;
		}

		std::string getTracebackAsSequence(const std::string& delim)
		{
			std::string seq = "";
			for (auto nodePtr : m_traceback_node_ptrs)
			{
				seq += nodePtr->getSequence() + delim;
			}
			return seq;
		}


		std::string getNormalizedCigarString()
		{
			return this->m_normalized_cigar_string;
		}

	private:

		void processTraceback(gssw_graph_mapping* graphMapping, size_t alignmentLength)
		{
			this->m_cigar_string = "";
			uint32_t totalScore = 0;
			m_traceback_node_ptrs.clear();
			m_soft_clip_occurrences = 0;
			uint32_t totalSoftclipLength = 0;
			gssw_node_cigar* nc = graphMapping->cigar.elements;
			for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
			{
				gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
				Node* nodePtr = (Node*)gsswNode->data;
				int32_t nodeScore = 0;
				uint32_t nodeLength = 0;
				uint32_t nodeSoftclipLength = 0;
				for (int j = 0; j < nc->cigar->length; ++j)
				{
					nodeLength += nc->cigar->elements[j].length;
					m_cigar_string += nc->cigar->elements[j].type + std::to_string(nc->cigar->elements[j].length);
					switch (nc->cigar->elements[j].type)
					{
					case 'M':
						nodeScore += (m_match_value * nc->cigar->elements[j].length);
						break;
					case 'X':
						nodeScore -= (m_mismatch_value * nc->cigar->elements[j].length);
						break;
					case 'I': // I and D are treated the same, this will fall through
					case 'D':
						nodeScore -= m_gap_open_value;
						nodeScore -= (m_gap_extension_value * (nc->cigar->elements[j].length -1));
						break;
					case 'S':
						nodeSoftclipLength += nc->cigar->elements[j].length;
						++m_soft_clip_occurrences;
					default:
						break;
					}
				}
				nodeScore = (nodeScore < 0) ? 0 : nodeScore; // the floor of the mapping score is 0
				totalSoftclipLength += nodeSoftclipLength;
				int32_t nodeScorePercent = (nodeLength > 0) ? ((float)nodeScore / ((float)(nodeLength - nodeSoftclipLength) * m_match_value)) * 100 : 0;
				totalScore += nodeScore;
				m_node_id_node_score_map.emplace(nodePtr->getID(), nodeScorePercent);
				m_traceback_node_ptrs.emplace_back(nodePtr);
			}
			this->m_total_score = ((float)totalScore / (float)((alignmentLength - totalSoftclipLength) * m_match_value)) * 100;
		}

		void setNormalizedCigarString(gssw_graph_mapping* graphMapping)
		{
			this->m_normalized_cigar_string = "";
			std::vector< char > cigTypes;
			std::vector< int > cigLens;
			gssw_node_cigar* nc = graphMapping->cigar.elements;
			for (auto i = 0; i < graphMapping->cigar.length; ++i, ++nc)
			{
				for (auto j = 0; j < nc->cigar->length; ++j)
				{
					if (cigTypes.size() > 0 && cigTypes[cigTypes.size() - 1] == nc->cigar->elements[j].type)
					{
						cigLens[cigTypes.size() - 1] += nc->cigar->elements[j].length;
					}
					else
					{
						cigTypes.emplace_back(nc->cigar->elements[j].type);
						cigLens.emplace_back(nc->cigar->elements[j].length);
					}
				}
			}
			for (auto i = 0; i < cigTypes.size(); ++i)
			{
				this->m_normalized_cigar_string += std::to_string(cigLens[i]) + cigTypes[i];
			}
		}

		std::string m_normalized_cigar_string;
		uint32_t m_soft_clip_occurrences;
		uint32_t m_total_score;
		std::string m_cigar_string;
		uint32_t m_match_value;
		uint32_t m_mismatch_value;
		uint32_t m_gap_open_value;
		uint32_t m_gap_extension_value;
		Graph::SharedPtr m_graph_ptr;
		std::vector< Node* > m_traceback_node_ptrs;
		std::unordered_map< uint32_t, uint32_t > m_node_id_node_score_map;
	};
}
