#ifndef GRAPHITE_VG_IVARIANTNODE
#define GRAPHITE_VG_IVARIANTNODE

#include "core/graph/INode.h"

#include "core/variant/Variant.h"
#include "core/util/Types.h"

namespace graphite
{
	namespace vg
	{
		class SNPNode;
		class IVariantNode : public INode
		{
		public:
			typedef std::shared_ptr< IVariantNode > SharedPtr;
		    IVariantNode(IVariant::SharedPtr variant, uint32_t altIndex) :
			INode(variant->getAltAllelePtrs()[altIndex]->getSequence(), strlen(variant->getAltAllelePtrs()[altIndex]->getSequence())),
				m_variant(variant),
				m_alt_index(altIndex)
			{
				m_position = variant->getPosition();
			}

			~IVariantNode() {}

			const char* getSequence() override { return m_variant->getAltAllelePtrs()[m_alt_index]->getSequence(); }

			static IVariantNode::SharedPtr BuildVariantNodes(IVariant::SharedPtr variant, uint32_t altIndex, const std::string& referenceSequence)
			{
				return std::dynamic_pointer_cast< IVariantNode >(std::make_shared< SNPNode >(variant, altIndex, referenceSequence));
			}

		private:
			IVariant::SharedPtr m_variant;
			size_t m_alt_index;
		};
	}// end namespace vg
}// end namespace graphite

#endif //GRAPHITE_VG_IVARIANTNODE
