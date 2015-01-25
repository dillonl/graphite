#ifndef GWIZ_VG_REFERENCE_NODE
#define GWIZ_VG_REFERENCE_NODE

#include "core/graph/INode.h"

namespace gwiz
{
	namespace vg
	{
		class ReferenceNode : public INode
		{
		public:
			typedef std::shared_ptr< ReferenceNode > SharedPtr;
		    ReferenceNode(IReference::SharedPtr referencePtr, char* sequence, uint32_t length) :
			    m_reference_ptr(referencePtr),
				INode(sequence, length)
				{}
			~ReferenceNode() {}

		private:
			IReference::SharedPtr m_reference_ptr;
		};
	}
}

#endif //GWIZ_VG_REFERENCE_NODE
