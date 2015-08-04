#include "GSSWMapping.h"


namespace gwiz
{
namespace gssw
{
    GSSWMapping::GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr) :
		m_gssw_mapping_ptr(gsswMappingPtr),
		m_alignment_ptr(alignmentPtr)
	{
		gssw_node_cigar* nc = m_gssw_mapping_ptr->cigar.elements;
		for (int i = 0; i < m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			m_allele_ptrs.emplace_back(((IAllele*)nc->node->data)->getSharedPtr());
		}
	}

	GSSWMapping::~GSSWMapping()
	{
	}

	int GSSWMapping::getMappingScore()
	{
		return this->m_gssw_mapping_ptr->score;
	}

	IAlignment::SharedPtr GSSWMapping::getAlignmentPtr()
	{
		return this->m_alignment_ptr;
	}

	std::vector< IAllele::SharedPtr > GSSWMapping::getAllelePtrs()
	{
		return this->m_allele_ptrs;
	}

}
}
