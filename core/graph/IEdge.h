#ifndef GRAPHITE_IEDGE_H
#define GRAPHITE_IEDGE_H

#include <boost/noncopyable.hpp>

#include <list>
#include <memory>

#include "INode.h"

namespace graphite
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

#endif // GRAPHITE_IEDGE_H
