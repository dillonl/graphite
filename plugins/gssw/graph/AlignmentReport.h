#ifndef GWIZ_GSSW_ALIGNMENTREPORT_H
#define GWIZ_GSSW_ALIGNMENTREPORT_H

#include "core/reference/IReference.h"
#include "core/variant/IVariantList.h"
#include "core/alignment/IAlignment.h"
#include "gssw/gssw.h"

#include <boost/noncopyable.hpp>
#include <memory>

namespace gwiz
{
namespace gssw
{
	class AlignmentReport : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< AlignmentReport > SharedPtr;

		AlignmentReport(IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignment::SharedPtr alignmentPtr, std::shared_ptr< gssw_graph_mapping > graphMappingPtr, position graphStartPosition);
		~AlignmentReport();

		std::string toString();

	private:
		IReference::SharedPtr m_reference_ptr;
		IVariantList::SharedPtr m_variant_list_ptr;
		IAlignment::SharedPtr m_alignment_ptr;
		std::shared_ptr< gssw_graph_mapping > m_graph_mapping_ptr;
		position m_graph_start_position;
	};
}
}

#endif //GWIZ_GSSW_ALIGNMENTREPORT_H
