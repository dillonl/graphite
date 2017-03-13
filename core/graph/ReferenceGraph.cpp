#include "ReferenceGraph.h"

namespace graphite
{
	ReferenceGraph::ReferenceGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, position startPosition, size_t graphSize, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue) :
		GSSWGraph(referencePtr, variantListPtr, startPosition, graphSize, matchValue, misMatchValue, gapOpenValue, gapExtensionValue)
	{
	}

	ReferenceGraph::~ReferenceGraph() {}

	void ReferenceGraph::constructGraph()
	{
		size_t graphStartOffset = this->m_start_position - this->m_reference_ptr->getRegion()->getStartPosition();

		std::string referenceSequenceString = std::string(this->m_reference_ptr->getSequence() + graphStartOffset, this->m_graph_size);

		auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
		auto referenceNodePtr = GSSWGraph::gssw_node_create_alt(this->m_start_position, referenceAllelePtr->getSequence(), referenceAllelePtr->getLength(), referenceAllelePtr, true, this->m_nt_table, this->m_mat);
		gssw_graph_add_node(this->m_graph_ptr, referenceNodePtr);
	}

	GSSWGraph::GSSWGraphMappingPtr ReferenceGraph::traceBackAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		gssw_graph_fill(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength(), this->m_nt_table, this->m_mat, this->m_gap_open, this->m_gap_extension, 15, 2);
		gssw_graph_mapping* graphMapping = gssw_graph_trace_back(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength(),m_match,m_mismatch,m_gap_open,m_gap_extension);
		gssw_node_cigar* nc = graphMapping->cigar.elements;
		for (int i = 0; i < graphMapping->cigar.length; ++i, ++nc)
		{
			nc->node->cigar = nc->cigar;
		}
		// gssw_print_graph_mapping(graphMapping);
		// gssw_graph_print_score_matrices(this->m_graph_ptr, alignmentPtr->getSequence(), alignmentPtr->getLength());
		auto graphMappingDeletor = [](gssw_graph_mapping* gm) { gssw_graph_mapping_destroy(gm); };
		return std::shared_ptr< gssw_graph_mapping >(graphMapping, graphMappingDeletor);
	}

}
