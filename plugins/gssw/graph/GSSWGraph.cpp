#include "GSSWGraph.h"
#include "vg/graph/ReferenceNode.h"

namespace gwiz
{
namespace gssw
{
	GSSWGraph::GSSWGraph(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr) :
		IGraph(referencePtr, variantListPtr), m_match(2), m_mismatch(2), m_gap_open(3), m_gap_extension(1)

	{
		this->m_nt_table = gssw_create_nt_table();
		this->m_mat = gssw_create_score_matrix(this->m_match, this->m_mismatch);
		this->m_graph_ptr = gssw_graph_create(0);
	}

	GSSWGraph::~GSSWGraph()
	{
		gssw_graph_destroy(this->m_graph_ptr);
		free(this->m_nt_table);
		free(this->m_mat);
	}

	void GSSWGraph::constructGraph()
	{
		position startPosition = this->m_reference_ptr->getRegion()->getStartPosition();
		size_t referenceOffset = 0;
		size_t referenceSize;
		Variant::SharedPtr variantPtr;
		std::vector< gssw_node* > altAndRefVertices;
		while (getNextCompoundVariant(variantPtr))
		{
			referenceSize = variantPtr->getPosition() - (startPosition + referenceOffset);
			if (referenceSize > 0)
			{
				auto referenceNode = gssw_node_create(this->m_reference_ptr->getSequence() + referenceOffset, referenceSize, this->m_nt_table, this->m_mat);
				gssw_graph_add_node(this->m_graph_ptr, referenceNode);
				altAndRefVertices.clear();
				altAndRefVertices.push_back(referenceNode);
			}

		}

		graphConstructed();
	}
}
}
