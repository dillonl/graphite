#ifndef GRAPHITE_GSSW_MAPPING_H
#define GRAPHITE_GSSW_MAPPING_H

#include "core/mapping/IMapping.h"

#include "gssw/gssw.h"

namespace graphite
{
namespace gssw
{

	class GSSWMapping : public IMapping
	{
    public:
		typedef std::shared_ptr< GSSWMapping > SharedPtr;
        GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr);
		~GSSWMapping();

		int getMappingScore();
		IAlignment::SharedPtr getAlignmentPtr();
		std::vector< IAllele::SharedPtr > getAllelePtrs();
		gssw_align* getGSSWAlignmentPtrFromAllelePtr(IAllele::SharedPtr allelePtr);

		void printLongFormat();

    private:

		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		std::unordered_map< IAllele::SharedPtr, gssw_align* > m_allele_alignment_map;
		IAlignment::SharedPtr m_alignment_ptr;
		uint32_t m_node_count;
	};

}
}

#endif //GRAPHITE_GSSW_MAPPING_H
