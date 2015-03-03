#ifndef GWIZ_GSSW_GSSWGRAPHMANAGER_H
#define GWIZ_GSSW_GSSWGRAPHMANAGER_H

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

		GraphManager(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, IAlignmentReader::SharedPtr alignmentReader, size_t padding);
		~GraphManager() {}

	private:
		void buildGraphs();
		std::queue< GSSWGraph::SharedPtr > m_gssw_graphs;

		gwiz::IReference::SharedPtr m_reference_ptr;
		gwiz::IVariantList::SharedPtr m_variant_list_ptr;
		IAlignmentReader::SharedPtr m_alignment_reader;
		size_t m_padding;
	};
}
}

#endif //GWIZ_GSSW_GSSWGRAPHMANAGER_H
