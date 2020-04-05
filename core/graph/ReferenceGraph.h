#ifndef GRAPHITE_REFERENCEGRAPH_H
#define GRAPHITE_REFERENCEGRAPH_H

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/alignment/Alignment.h"

#include "Node.h"

#include "gssw.h"

#include <memory>

namespace graphite
{

	class ReferenceGraph : private Noncopyable
	{
	public:
		typedef std::shared_ptr< ReferenceGraph > SharedPtr;
		ReferenceGraph(const std::string& refSequence, position startPosition);
		~ReferenceGraph();

		int32_t adjudicateAlignment(Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t gapExtensionValue);

	private:
		int32_t processTraceback(gssw_graph_mapping* graphMapping, Alignment::SharedPtr alignmentPtr, Sample::SharedPtr samplePtr, bool isForwardStrand, uint32_t  matchValue, uint32_t mismatchValue, uint32_t gapOpenValue, uint32_t  gapExtensionValue);

		Node::SharedPtr m_node;
		Region::SharedPtr m_region_ptr;
		std::unordered_set< std::string > m_aligned_read_names;
		std::mutex m_aligned_read_names_mutex;

	};

}

#endif //GRAPHITE_REFERENCEGRAPH_H
