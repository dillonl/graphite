#ifndef GWIZ_VG_IVARIANTNODE
#define GWIZ_VG_IVARIANTNODE

#include "core/graph/INode.h"
#include "core/variants/Variant.h"
#include "core/utils/Types.h"

namespace gwiz
{
	namespace vg
	{
		class IVariantNode : public INode
		{
		public:
			typedef std::shared_ptr< IVariantNode > SharedPtr;
		    IVariantNode(IVariant::SharedPtr variant) : m_variant(variant)
			{}
			IVariantNode() {}
			~IVariantNode() {}

			IVariant::SharedPtr m_variant;

			std::vector< IVariantNode::SharedPtr > BuildVariantNodes(Variant::SharedPtr variant)
			{
				std::vector< IVariantNode::SharedPtr > variantNodes;
				for (auto iter = variant->getAlt().begin(), iter != variant->getAlt().end(); ++iter)
				{
					/* variantNodes.push_back(std::make_shared< SNPNode >()); */
				}
				return variantNodes;
			}
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_IVARIANTNODE
