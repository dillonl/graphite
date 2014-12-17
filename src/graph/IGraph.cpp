#include "IGraph.h"
#include "graph/SparseGraph.h"

namespace gwiz
{
	std::shared_ptr<IGraph> IGraph::BuildGraph(std::shared_ptr<IReference> reference, std::list< std::shared_ptr<IVariant> >, GRAPH_TYPE graphType)
		{
			std::shared_ptr<IGraph> graph(NULL);
			switch (graphType)
			{
			case GRAPH_TYPE::SPARSE:
				graph = std::make_shared<SparseGraph>();
				break;
			case GRAPH_TYPE::DENSE:
				//graph = std::make_shared<DenseGraph>();
				break;
			default:
				break;
			}
			graph->constructGraph();

			return graph;
		}
}
