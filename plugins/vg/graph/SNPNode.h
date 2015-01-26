#ifndef GWIZ_VG_SNPNODE
#define GWIZ_VG_SNPNODE

#include "core/graph/INode.h"
#include "core/utils/Types.h"

namespace gwiz
{
	namespace vg
	{
		class SNPNode : public INode
		{
		public:
			typedef std::shared_ptr< SNPNode > SharedPtr;
		    SNPNode(const char* sequence, position position, uint32_t length) :
			    INode(sequence, position, length)
			    {}
			~SNPNode() {}

			std::string m_test;
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_SNPNODE
