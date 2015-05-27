#ifndef GWIZ_VG_SNPNODE
#define GWIZ_VG_SNPNODE

#include "IVariantNode.h"
#include "core/graph/INode.h"
#include "core/util/Types.h"

namespace gwiz
{
	namespace vg
	{
		class SNPNode : public IVariantNode
		{
		public:
			typedef std::shared_ptr< SNPNode > SharedPtr;

    		SNPNode(Variant::SharedPtr variant, uint32_t altIndex) :
				IVariantNode(variant, altIndex)
			{
			}
			~SNPNode() {}

		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_SNPNODE
