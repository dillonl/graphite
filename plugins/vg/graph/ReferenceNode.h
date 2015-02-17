#ifndef GWIZ_VG_REFERENCE_NODE
#define GWIZ_VG_REFERENCE_NODE

#include "core/reference/IReference.h"
#include "core/graph/INode.h"

namespace gwiz
{
	namespace vg
	{
		class ReferenceNode : public INode
		{
		public:
			typedef std::shared_ptr< ReferenceNode > SharedPtr;
    		ReferenceNode(IReference::SharedPtr reference, size_t offset, size_t length) :
				m_offset(offset),
				INode(reference->getSequence() + offset, length)
				{
					m_position = reference->getRegion()->getStartPosition() + offset;
				}
			~ReferenceNode() {}



		private:
			size_t m_offset;
		};
	}
}

#endif //GWIZ_VG_REFERENCE_NODE
