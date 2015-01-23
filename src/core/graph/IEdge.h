#ifndef GWIZ_IEDGE_H
#define GWIZ_IEDGE_H

#include <list>
#include <memory>

#include "core/utils/Noncopyable.h"
#include "INode.h"

namespace gwiz
{
	class IEdge : noncopyable
	{
	public:
		typedef std::shared_ptr<IEdge> SharedPtr;
		IEdge() {}
		virtual ~IEdge() {}

	protected:

	private:

	};
}

#endif // GWIZ_IEDGE_H
