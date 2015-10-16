#ifndef GRAPHITE_VG_SNPNODE
#define GRAPHITE_VG_SNPNODE

#include "IVariantNode.h"
#include "core/graph/INode.h"
#include "core/util/Types.h"

namespace graphite
{
	namespace vg
	{
		class SNPNode : public IVariantNode
		{
		public:
			typedef std::shared_ptr< SNPNode > SharedPtr;

    		SNPNode(IVariant::SharedPtr variantPtr, uint32_t altIndex) :
				IVariantNode(variantPtr, altIndex)
			{
			}
			~SNPNode() {}

		};
	}// end namespace vg
}// end namespace graphite

#endif //GRAPHITE_VG_SNPNODE
