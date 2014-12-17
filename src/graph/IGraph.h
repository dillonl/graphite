#ifndef GWIZ_IGRAPH_H
#define GWIZ_IGRAPH_H

#include "core/IReference.h"
#include "core/IVariant.h"
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
		IGraph() {}
		virtual ~IGraph() {}


		static std::shared_ptr<IGraph> BuildGraph(std::shared_ptr<IReference> reference, std::list< std::shared_ptr<IVariant> >, GRAPH_TYPE graphType);
	protected:

		virtual void constructGraph() = 0;

	private:


	};
} // end namespace gwiz

#endif // GWIZ_IGRAPH_H
