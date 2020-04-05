#include "ReferenceGraph.h"


namespace graphite
{

	ReferenceGraph::ReferenceGraph(const std::string& refSequence, position startPosition)
	{
		this->m_node = std::make_shared< Node >(refSequence, startPosition, Node::ALLELE_TYPE::REF);
	}

	ReferenceGraph::~ReferenceGraph()
	{
	}

	int32_t ReferenceGraph::adjudicateAlignment(Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue)
	{
		gssw_sse2_disable();
		int8_t* nt_table = gssw_create_nt_table();
		int8_t* mat = gssw_create_score_matrix(matchValue, mismatchValue);
		gssw_graph* graph = gssw_graph_create(1);

		gssw_node* gsswNode = (gssw_node*)gssw_node_create(m_node.get(), m_node->getID(), m_node->getSequence().c_str(), nt_table, mat);
		gssw_graph_add_node(graph, gsswNode);

		gssw_graph_fill(graph, alignmentPtr->getSequence(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0, 15, 2, true);
		gssw_graph_mapping* gm = gssw_graph_trace_back (graph, alignmentPtr->getSequence(), alignmentPtr->getLength(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0);
		float swPercent = processTraceback(gm, alignmentPtr, samplePtr, alignmentPtr->getIsForwardStrand(), matchValue, mismatchValue, gapOpenValue, gapExtensionValue);
		gssw_graph_mapping_destroy(gm);

		// note that nodes which are referred to in this graph are destroyed as well
		gssw_graph_destroy(graph);

		free(nt_table);
		free(mat);
		return swPercent;
	}

	int32_t  ReferenceGraph::processTraceback(gssw_graph_mapping* graphMapping, Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue)
	{
		int count = 0;

		std::vector< std::tuple< Node*, uint32_t > > nodePtrScoreTuples;
		uint32_t totalScore = 0;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		uint32_t softclipLength = 0;
		std::string fullCigarString = "";
		uint32_t softclipCount = 0;
		bool hasAlternate = false;
		uint32_t totalSoftclipLength = 0;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
			Node* nodePtr = (Node*)gsswNode->data;
			int32_t nodeScore = 0;
			uint32_t nodeLength = 0;
			uint32_t nodeSoftclipLength = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				std::tuple< uint32_t, char > cigarComponent = std::make_tuple(nc->cigar->elements[j].length, nc->cigar->elements[j].type);
				nodeLength += nc->cigar->elements[j].length;
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					nodeScore += (matchValue * nc->cigar->elements[j].length);
					break;
				case 'X':
					nodeScore -= (mismatchValue * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					nodeScore -= gapOpenValue;
					nodeScore -= (gapExtensionValue * (nc->cigar->elements[j].length -1));
					break;
				case 'S':
					nodeSoftclipLength += nc->cigar->elements[j].length;
				default:
					break;
				}
			}
			nodeScore = (nodeScore < 0) ? 0 : nodeScore; // the floor of the mapping score is 0
			totalSoftclipLength += nodeSoftclipLength;
			int32_t nodeScorePercent = (nodeLength > 0) ? ((float)nodeScore / ((float)(nodeLength - nodeSoftclipLength) * matchValue)) * 100 : 0;
			totalScore += nodeScore;
		}

		int32_t totalScorePercent = ((float)(totalScore))/((float)(alignmentPtr->getLength() - softclipLength)) * 100;
		return totalScorePercent;
	}

}
