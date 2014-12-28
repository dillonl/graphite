#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/reference/IReference.h"
#include "core/variants/IVariant.h"
#include "utils/NonCopyable.h"

#include <list>
#include <string>
#include <memory>

namespace gwiz
{
	enum class GRAPH_TYPE {SPARSE, DENSE};

	class IGraph : private noncopyable
	{
	public:
		typedef std::shared_ptr<IGraph> SharedPtr;

		IGraph() {}
		virtual ~IGraph() {}


		static IGraph::SharedPtr BuildGraph(IReference::SharedPtr reference, std::list< IVariant::SharedPtr >, GRAPH_TYPE graphType);
	protected:

		virtual void constructGraph() = 0;

	private:


	};
} // end namespace gwiz

#endif // GWIZ_IGRAPH_H
