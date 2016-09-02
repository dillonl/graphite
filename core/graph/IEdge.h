#ifndef GRAPHITE_IEDGE_H
#define GRAPHITE_IEDGE_H

#include <list>
#include <memory>

#include "INode.h"
#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class IEdge : Noncopyable
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
