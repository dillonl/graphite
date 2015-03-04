#ifndef GWIZ_GSSW_GSSWGRAPHMANAGER_H
#define GWIZ_GSSW_GSSWGRAPHMANAGER_H

#include "core/alignments/IAlignmentReaderManager.h"
#include "GSSWGraph.h"

#include <queue>
#include <memory>

#include <boost/noncopyable.hpp>

namespace gwiz
{
namespace gssw
{
	class GraphManager : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< GraphManager > SharedPtr;

		GraphManager(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, size_t padding);
		~GraphManager() {}

		void buildGraphs();

	private:
		std::queue< GSSWGraph::SharedPtr > m_gssw_graphs;

		gwiz::IReference::SharedPtr m_reference_ptr;
		gwiz::IVariantList::SharedPtr m_variant_list_ptr;
		IAlignmentReaderManager::SharedPtr m_alignment_reader_manager;
		size_t m_padding;
	};
}
}

#endif //GWIZ_GSSW_GSSWGRAPHMANAGER_H
