#include "VariantContig.h"

#include <map>

namespace gwiz
{
namespace adjudicator
{

	VariantContig::VariantContig(Region::SharedPtr region, IReference::SharedPtr referencePtr, vg::ReferenceNode::SharedPtr startReferenceNode, vg::ReferenceNode::SharedPtr endReferenceNode) :
		m_region(region),
		m_reference_ptr(referencePtr),
		m_start_reference_node(startReferenceNode),
		m_end_reference_node(endReferenceNode)
	{
	}

	VariantContig::~VariantContig()
	{
	}

	void VariantContig::buildVariantContig()
	{
		this->m_contigs.clear();

	}

} // namespace adjudicator
} // namespace gwiz
