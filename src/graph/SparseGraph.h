#ifndef GWIZ_SPARSE_GRAPH_H
#define GWIZ_SPARSE_GRAPH_H

#include "graph/IGraph.h"

namespace gwiz
{
	class SparseGraph : public IGraph
	{
	public:
		SparseGraph();
		~SparseGraph();

	protected:
		void constructGraph() override;
	private:

	};
}

#endif // GWIZ_SPARSE_GRAPH_H
