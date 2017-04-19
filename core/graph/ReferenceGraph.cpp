#include "ReferenceGraph.h"

namespace graphite
{
	ReferenceGraph::ReferenceGraph(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, Region::SharedPtr regionPtr, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue, uint32_t numGraphCopies) :
		GSSWGraph(referencePtr, variantListPtr, regionPtr, matchValue, misMatchValue, gapOpenValue, gapExtensionValue, numGraphCopies)
	{
	}

	ReferenceGraph::~ReferenceGraph() {}

	void ReferenceGraph::constructGraph()
	{
		std::string referenceSequenceString = this->m_reference_ptr->getSequenceFromRegion(this->m_region_ptr);

		auto referenceAllelePtr = std::make_shared< Allele >(referenceSequenceString);
		auto referenceNodePtr = GSSWGraph::gssw_node_create_alt(this->m_region_ptr->getStartPosition(), referenceAllelePtr->getSequence(), referenceAllelePtr->getLength(), referenceAllelePtr, true, this->m_nt_table, this->m_mat);
		gssw_graph_add_node(this->m_graph_ptr, referenceNodePtr);
		generateGraphCopies();
	}

}
