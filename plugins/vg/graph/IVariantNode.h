#ifndef GWIZ_VG_IVARIANTNODE
#define GWIZ_VG_IVARIANTNODE

#include "core/graph/INode.h"

#include "core/variant/Variant.h"
#include "core/util/Types.h"

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
			INode(variant->getAltAllelePtrs()[altIndex]->getSequence(), strlen(variant->getAltAllelePtrs()[altIndex]->getSequence())),
				m_variant(variant),
				m_alt_index(altIndex)
			{
				m_position = variant->getPosition();
			}

			~IVariantNode() {}

			const char* getSequence() override { return m_variant->getAltAllelePtrs()[m_alt_index]->getSequence(); }

			static IVariantNode::SharedPtr BuildVariantNodes(Variant::SharedPtr variant, uint32_t altIndex)
			{
				return std::dynamic_pointer_cast< IVariantNode >(std::make_shared< SNPNode >(variant, altIndex));
			}

		private:
			Variant::SharedPtr m_variant;
			size_t m_alt_index;
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_IVARIANTNODE
