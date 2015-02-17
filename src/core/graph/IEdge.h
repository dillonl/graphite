#ifndef GWIZ_IEDGE_H
#define GWIZ_IEDGE_H

#include <boost/noncopyable.hpp>

#include <list>
#include <memory>

#include "INode.h"

namespace gwiz
{
	class IEdge : boost::noncopyable
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
