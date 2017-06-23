#ifndef GRAPHITE_ADJUDICATOR_MAPPING_H
#define GRAPHITE_ADJUDICATOR_MAPPING_H

#include "core/mapping/IMapping.h"
#include "core/adjudicator/IAdjudicator.h"
#include "core/graph/GraphManager.h"

#include "gssw.h"

namespace graphite
{

	class GSSWMapping : public IMapping
	{
    public:
		typedef std::shared_ptr< GSSWMapping > SharedPtr;
        GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr);
		~GSSWMapping();

		int getMappingScore() override;
		/* MappingAlignmentInfo::SharedPtr getMappingAlignmentInfo(IAllele::SharedPtr allelePtr, IAdjudicator::SharedPtr adjudicatorPtr) override; */
		IAlignment::SharedPtr getAlignmentPtr() override;
		std::vector< IAllele::SharedPtr > getAllelePtrs() override;
		position getPosition() override { return m_position; }
        uint32_t getOffset () 
        { 
            if (m_offset == 0)
            {
                m_offset = 1;
                return m_offset;
            }

            // Offset from the beginning of the graph. +1 because the started base will be at position 1.
            return m_offset + 1;
        }      
		std::vector< MappingAlignmentInfo::SharedPtr > getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr);
		void incrementAlleleCounts() override;
		void setMapped(bool mapped) override;
		bool getMapped() override { return m_mapped; }
		void addAlleleCountCallback(std::function< void () > functor) override;

		void printMapping() override;
		void printSimpleMapping();

        void getGraphPathHeaderAndSequence (std::string& graphPathHeader, std::string& graphPathSequence, position variantPosition);  // Append the graph path headers to each node during traceback.

        std::string getCigarString (IAdjudicator::SharedPtr adjudicatorPtr);

        std::vector< uint32_t > getNodeIDs ();
        std::vector< int32_t > getNodeLengths ();
        std::vector< uint32_t > getNodeStartPositions ();
        std::vector< uint32_t > getNodeEndPositions ();
        std::vector< std::string > getRefOrAltStrings ();
        std::vector< NodeInfo::VariantType > getVariantTypes ();

    private:

		std::shared_ptr< gssw_graph_mapping > m_gssw_mapping_ptr;
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		std::unordered_map< IAllele*, gssw_node* > m_allele_gssw_nodes_map;
		IAlignment::SharedPtr m_alignment_ptr;
		position m_position;
		std::vector< std::function< void () > > m_allele_incrementor_callback_list;
		bool m_mapped;

        uint32_t m_offset;      // Offset from the starting of the graph.
        
        // BED info.
        std::string m_graph_path_header;
        std::vector< uint32_t > m_node_IDs;
        std::vector< int32_t > m_node_lengths;
        /*
        std::vector< uint32_t > m_node_start_positions;
        std::vector< uint32_t > m_node_end_positions;
        std::vector< std::string > m_ref_or_alt_strings;
        */
        std::vector< NodeInfo::VariantType > m_variant_types;

	};

}

#endif //GRAPHITE_GSSW_MAPPING_H
