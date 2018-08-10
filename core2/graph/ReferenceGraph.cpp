#include "ReferenceGraph.h"


namespace graphite
{
	ReferenceGraph::ReferenceGraph(FastaReference::SharedPtr fastaReferencePtr, Region::SharedPtr regionPtr) :
		m_region_ptr(regionPtr)
	{
		std::string referenceSequence = fastaReferencePtr->getSequenceStringFromRegion(regionPtr);
		this->m_node = std::make_shared< Node >(referenceSequence, regionPtr->getStartPosition(), Node::ALLELE_TYPE::REF);
	}

	ReferenceGraph::~ReferenceGraph()
	{
	}

	float ReferenceGraph::adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue)
	{
		gssw_sse2_disable();
		int8_t* nt_table = gssw_create_nt_table();
		int8_t* mat = gssw_create_score_matrix(matchValue, mismatchValue);
		gssw_graph* graph = gssw_graph_create(1);

		gssw_node* gsswNode = (gssw_node*)gssw_node_create(m_node.get(), m_node->getID(), m_node->getSequence().c_str(), nt_table, mat);
		gssw_graph_add_node(graph, gsswNode);

		gssw_graph_fill(graph, bamAlignmentPtr->QueryBases.c_str(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0, 15, 2, true);
		gssw_graph_mapping* gm = gssw_graph_trace_back (graph, bamAlignmentPtr->QueryBases.c_str(), bamAlignmentPtr->QueryBases.size(), nt_table, mat, gapOpenValue, gapExtensionValue, 0, 0);
		float swPercent = processTraceback(gm, bamAlignmentPtr, samplePtr, !bamAlignmentPtr->IsReverseStrand(), matchValue, mismatchValue, gapOpenValue, gapExtensionValue);
		gssw_graph_mapping_destroy(gm);

		// note that nodes which are referred to in this graph are destroyed as well
		gssw_graph_destroy(graph);

		free(nt_table);
		free(mat);
		return swPercent;
	}

	float ReferenceGraph::processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue)
	{
		std::string alignmentName = bamAlignmentPtr->Name + std::to_string(bamAlignmentPtr->IsFirstMate());
		{
			std::lock_guard< std::mutex > l(m_aligned_read_names_mutex);
			if (this->m_aligned_read_names.find(alignmentName) != this->m_aligned_read_names.end())
			{
				return 0;
			}
			this->m_aligned_read_names.emplace(alignmentName);
		}
		std::vector< std::tuple< Node*, uint32_t > > nodePtrScoreTuples;
		uint32_t prefixMatch = 0;
		uint32_t suffixMatch = 0;
		bool setPrefix = true;
		uint32_t totalScore = 0;
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		uint32_t softclipLength = 0;
		std::string fullCigarString = "";
		uint32_t softclipCount = 0;
		bool hasAlternate = false;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			std::string cigarString = "";
			gssw_node* gsswNode = graphMapping->cigar.elements[i].node;
			Node* nodePtr = (Node*)gsswNode->data;
			hasAlternate |= (nodePtr->getAlleleType() == Node::ALLELE_TYPE::ALT);
			int32_t score = 0;
			uint32_t length = 0;
			uint32_t tmpSofclipLength = 0;
			// std::unordered_map< gssw_node*, uint32_t > nodePrefixMatchMap;
			// std::unordered_map< gssw_node*, uint32_t > nodeSuffixMatchMap;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				// std::cout << nc->cigar->elements[j].length << nc->cigar->elements[j].type;
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					suffixMatch += nc->cigar->elements[j].length;
					if (setPrefix)
					{
						prefixMatch += nc->cigar->elements[j].length;
					}
					score += (matchValue * nc->cigar->elements[j].length);
					break;
				case 'X':
					setPrefix = false;
					suffixMatch = 0;
					score -= (mismatchValue * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					setPrefix = false;
					suffixMatch = 0;
					score -= gapOpenValue;
					score -= (gapExtensionValue * (nc->cigar->elements[j].length -1));
					break;
				case 'S':
					tmpSofclipLength += nc->cigar->elements[j].length;
					++softclipCount;
				default:
					break;
				}
				length += nc->cigar->elements[j].length;
				fullCigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			// std::cout << std::endl;
			score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
			softclipLength += tmpSofclipLength;
			float tmpScore = score;
			uint32_t alleleSWScorePercent = (length > 0) ? (tmpScore / (length - tmpSofclipLength)) * 100 : 0;
			std::tuple< Node*, uint32_t > nodePtrScoreTuple(nodePtr, alleleSWScorePercent);
			nodePtrScoreTuples.emplace_back(nodePtrScoreTuple);
			totalScore += score;
		}

		float totalScorePercent = ((float)(totalScore))/((float)(bamAlignmentPtr->Length - softclipLength)) * 100;
		return totalScorePercent;
	}

}
