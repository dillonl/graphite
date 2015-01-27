#ifndef GWIZ_VG_IVARIANTNODE
#define GWIZ_VG_IVARIANTNODE

#include "core/graph/INode.h"

#include "core/variants/Variant.h"
#include "core/utils/Types.h"

namespace gwiz
{
	namespace vg
	{
		class SNPNode;
		class IVariantNode : public INode
		{
		public:
			typedef std::shared_ptr< IVariantNode > SharedPtr;
		    IVariantNode(Variant::SharedPtr variant, uint32_t altIndex) :
			    INode(variant->getAlt()[altIndex].c_str(), variant->getAlt()[altIndex].size()),
			    m_variant(variant)
			{}

			~IVariantNode() {}

			static IVariantNode::SharedPtr BuildVariantNodes(Variant::SharedPtr variant, uint32_t altIndex)
			{
				return std::dynamic_pointer_cast< IVariantNode >(std::make_shared< SNPNode >(variant, altIndex));
			}

		private:
			Variant::SharedPtr m_variant;
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_IVARIANTNODE
