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
		    ReferenceNode(const char* sequence, position position, uint32_t length) :
			    INode(sequence, position, length)
				{}
			~ReferenceNode() {}

			void updateLength(uint32_t length) { m_length = length; }
		};
	}
}

#endif //GWIZ_VG_REFERENCE_NODE
