#ifndef GWIZ_VG_VGNODE
#define GWIZ_VG_VGNODE

#include "core/graph/INode.h"
#include "core/util/Types.h"

namespace gwiz
{
	namespace vg
	{
		class VGNode : INode
		{
		public:
			typedef std::shared_ptr< VGNode > SharedPtr;
			VGNode();
			~VGNode();



		};
	}// end namespace vg
}// end namespace gwiz

#endif //GWIZ_VG_VGNODE
