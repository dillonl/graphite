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
		SNPNode(std::string& test) : m_test(test)
			{}
			~SNPNode() {}

			std::string m_test;
		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_SNPNODE
