#ifndef GWIZ_GSSW_MAPPING_H
#define GWIZ_GSSW_MAPPING_H

#include "core/mapping/IMapping.h"

#include "gssw/gssw.h"

namespace gwiz
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

    private:
		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		IAlignment::SharedPtr m_alignment_ptr;
	};

}
}

#endif //GWIZ_GSSW_MAPPING_H
