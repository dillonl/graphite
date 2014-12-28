#ifndef GWIZ_IEDGE_H
#define GWIZ_IEDGE_H

#include <list>
#include <memory>

#include "utils/Noncopyable.h"
#include "INode.h"

namespace gwiz
{
	class IEdge : noncopyable
	{
	public:
		typedef std::shared_ptr<IEdge> SharedPtr;

		virtual ~IEdge() {}

	protected:
		std::list<INode::SharedPtr> m_nodes;

	private:

	};
}

#endif // GWIZ_IEDGE_H
