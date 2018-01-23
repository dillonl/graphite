#ifndef GRAPHITE_ADJUDICATOR_MAPPING_H
#define GRAPHITE_ADJUDICATOR_MAPPING_H

#include "core/mapping/IMapping.h"
#include "core/adjudicator/IAdjudicator.h"

#include "gssw.h"

namespace graphite
{

	class GSSWMapping : public IMapping
	{
    public:
		typedef std::shared_ptr< GSSWMapping > SharedPtr;
        GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr, position startPosition);
		~GSSWMapping();

		int getMappingScore() override;
		IAlignment::SharedPtr getAlignmentPtr() override;
		std::vector< IAllele::SharedPtr > getAllelePtrs() override;
		position getPosition() override { return m_position; }
		std::vector< MappingAlignmentInfo::SharedPtr > getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr);
		void incrementAlleleCounts() override;
		void setMapped(bool mapped) override;
		bool getMapped() override { return m_mapped; }
		void addAlleleCountCallback(std::function< void () > functor) override;
		void setTracebackSequenceAndID(std::string& sequence, std::string& id) override;
		std::vector< std::tuple< char, uint32_t > > getCigarData() override;
		position getAlignmentMappedPosition() override { return m_offset; }

		void printMapping() override;
		void printSimpleMapping();

    private:

		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		std::unordered_map< IAllele*, gssw_node* > m_allele_gssw_nodes_map;
		IAlignment::SharedPtr m_alignment_ptr;
		position m_position;
		std::vector< std::function< void () > > m_allele_incrementor_callback_list;
		bool m_mapped;
		position m_offset;
	};

}

#endif //GRAPHITE_GSSW_MAPPING_H
