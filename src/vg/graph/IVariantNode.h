#ifndef GWIZ_VG_IVARIANTNODE
#define GWIZ_VG_IVARIANTNODE

#include "core/graph/INode.h"
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
			~IVariantNode() {}

			IVariant::SharedPtr m_variant;
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_IVARIANTNODE
