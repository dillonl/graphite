#ifndef GRAPHITE_ADJUDICATOR_ALIGNMENTREPORT_H
#define GRAPHITE_ADJUDICATOR_ALIGNMENTREPORT_H

#include "core/reference/IReference.h"
#include "core/variant/VariantList.h"
#include "core/alignment/IAlignment.h"
#include "core/util/Noncopyable.hpp"

#include "gssw.h"

#include <memory>

namespace graphite
{
	class AlignmentReport : private Noncopyable
	{
	public:
		typedef std::shared_ptr< AlignmentReport > SharedPtr;

		AlignmentReport(IReference::SharedPtr referencePtr, VariantList::SharedPtr variantListPtr, IAlignment::SharedPtr alignmentPtr, std::shared_ptr< gssw_graph_mapping > graphMappingPtr, position graphStartPosition);
		~AlignmentReport();

		std::string toString();

	private:
		IReference::SharedPtr m_reference_ptr;
		VariantList::SharedPtr m_variant_list_ptr;
		IAlignment::SharedPtr m_alignment_ptr;
		std::shared_ptr< gssw_graph_mapping > m_graph_mapping_ptr;
		position m_graph_start_position;
	};
}

#endif //GRAPHITE_GSSW_ALIGNMENTREPORT_H
