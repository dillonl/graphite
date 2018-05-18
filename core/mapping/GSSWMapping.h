#ifndef GRAPHITE_ADJUDICATOR_MAPPING_H
#define GRAPHITE_ADJUDICATOR_MAPPING_H

#include "core/adjudicator/GSSWAdjudicator.h"
#include "MappingAlignmentInfo.h"
#include "core/alignment/IAlignment.h"
#include "core/util/Noncopyable.hpp"

#include "gssw.h"

#include <memory>

namespace graphite
{

    class GSSWAdjudicator;
	class GSSWMapping : private Noncopyable
	{
    public:
		typedef std::shared_ptr< GSSWMapping > SharedPtr;
        GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr, position startPosition);
		~GSSWMapping();

		int getMappingScore();
		IAlignment::SharedPtr getAlignmentPtr();
		std::vector< IAllele::SharedPtr > getAllelePtrs();
		position getPosition() { return m_position; }
        std::vector< MappingAlignmentInfo::SharedPtr > getMappingAlignmentInfoPtrs(std::shared_ptr< GSSWAdjudicator > adjudicatorPtr);
		void incrementAlleleCounts();
		void setMapped(bool mapped);
		bool getMapped() { return m_mapped; }
		void addAlleleCountCallback(std::function< void () > functor);
		void setTracebackSequenceAndID(std::string& sequence, std::string& id);
		std::vector< std::tuple< char, uint32_t > > getCigarData();
		position getAlignmentMappedPosition() { return m_offset; }
		void printMapping();
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
