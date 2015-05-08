#ifndef GWIZ_GSSW_GSSWADJUDICATOR_H
#define GWIZ_GSSW_GSSWADJUDICATOR_H

#include "core/adjudicator/IGraphAdjudicator.h"
#include "GSSWGraph.h"

#include <mutex>

namespace gwiz
{
namespace gssw
{
	class GSSWAdjudicator : public IGraphAdjudicator
	{
	public:
		typedef std::shared_ptr< GSSWAdjudicator > SharedPtr;
		GSSWAdjudicator();
		~GSSWAdjudicator();

		IVariantList::SharedPtr adjudicateGraph(IGraph::SharedPtr graphPtr, IAlignmentReader::SharedPtr alignmentsPtr) override;
		int m_id;
	private:
		void printNodes(GSSWGraph::SharedPtr graphPtr, const std::string& alignment);
		bool isAmbiguousFinalNodeCall(gssw_node* node)
		{
			if (node->count_prev > 0)
			{
				auto prevNode = node->prev[0];
				gssw_node* tmpNode;
				/* for (int i = 0; i <  */
			}
			return false;
		}
		std::mutex m_adjudication_lock;
	};
}
}

#endif //GWIZ_GSSW_GSSWADJUDICATOR_H
