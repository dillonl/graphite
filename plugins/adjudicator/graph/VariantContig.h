#ifndef GWIZ_ADJUDICATOR_VARIANTCONTIG_H
#define GWIZ_ADJUDICATOR_VARIANTCONTIG_H

#include <boost/noncopyable.hpp>


#include "core/reference/IReference.h"
#include "vg/graph/ReferenceNode.h"

namespace gwiz
{
namespace adjudicator
{

	class VariantContig : boost::noncopyable
	{
	public:
        VariantContig(Region::SharedPtr region, IReference::SharedPtr referencePtr, vg::ReferenceNode::SharedPtr startReferenceNode, vg::ReferenceNode::SharedPtr endReferenceNode);
		~VariantContig();

		void buildVariantContig();

	private:
		Region::SharedPtr m_region;
		IReference::SharedPtr m_reference_ptr;
        vg::ReferenceNode::SharedPtr m_start_reference_node;
        vg::ReferenceNode::SharedPtr m_end_reference_node;
		std::vector< std::string > m_contigs;
	};

} // namespace adjudicator
} // namespace gwiz

#endif //GWIZ_ADJUDICATOR_VARIANTCONTIG_H
