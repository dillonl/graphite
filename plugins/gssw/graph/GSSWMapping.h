#ifndef GRAPHITE_GSSW_MAPPING_H
#define GRAPHITE_GSSW_MAPPING_H

#include "core/mapping/IMapping.h"
#include "core/adjudicator/IAdjudicator.h"

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
		MappingAlignmentInfo::SharedPtr getMappingAlignmentInfo(IAllele::SharedPtr allelePtr, IAdjudicator::SharedPtr adjudicatorPtr);
		IAlignment::SharedPtr getAlignmentPtr();
		std::vector< IAllele::SharedPtr > getAllelePtrs();
		position getPosition() override { return m_position; }
		std::vector< MappingAlignmentInfo::SharedPtr > getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr);
		void incrementAlleleCounts() override;
		void setMapped(bool mapped) override;
		bool getMapped() override { return m_mapped; }
		void addAlleleCountCallback(std::function< void () > functor);

		void printLongFormat();

    private:

		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		std::unordered_map< IAllele*, gssw_node* > m_allele_gssw_nodes_map;
		IAlignment::SharedPtr m_alignment_ptr;
		position m_position;
		std::vector< std::function< void () > > m_allele_incrementor_callback_list;
		bool m_mapped;
	};

}
}

#endif //GRAPHITE_GSSW_MAPPING_H
