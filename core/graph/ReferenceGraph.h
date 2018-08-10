#ifndef GRAPHITE_REFERENCEGRAPH_H
#define GRAPHITE_REFERENCEGRAPH_H

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"

#include "Node.h"

#include "api/BamAlignment.h"

#include "gssw.h"

#include <memory>

namespace graphite
{

	class ReferenceGraph : private Noncopyable
	{
	public:
		typedef std::shared_ptr< ReferenceGraph > SharedPtr;
		ReferenceGraph(FastaReference::SharedPtr fastaReferencePtr, Region::SharedPtr regionPtr);
		~ReferenceGraph();

float adjudicateAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue);

	private:
		float processTraceback(gssw_graph_mapping* graphMapping, std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue);

		Node::SharedPtr m_node;
		Region::SharedPtr m_region_ptr;
		std::unordered_set< std::string > m_aligned_read_names;
		std::mutex m_aligned_read_names_mutex;

	};

}

#endif //GRAPHITE_REFERENCEGRAPH_H
