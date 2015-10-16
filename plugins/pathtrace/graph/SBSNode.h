#ifndef GRAPHITE_VG_VGNODE
#define GRAPHITE_VG_VGNODE

#include "core/graph/INode.h"
#include "core/util/Types.h"

namespace graphite
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
}// end namespace graphite

#endif //GRAPHITE_VG_VGNODE
