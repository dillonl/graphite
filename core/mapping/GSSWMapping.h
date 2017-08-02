#ifndef GRAPHITE_ADJUDICATOR_MAPPING_H
#define GRAPHITE_ADJUDICATOR_MAPPING_H

#include "core/adjudicator/IAdjudicator.h"
#include "core/mapping/IMapping.h"
#include "core/mapping/NodeInfo.h"
//#include "core/graph/GraphManager.h"    // There is a circular dependency between GraphManager, SAMFileWriter and GSSWMapping. I'm pretty sure the only thing that GSSWMapping needs from GraphManager is the NodeInfo class. Find out the best way to address this problem when I get the SAMFileWriter where I want it in the GraphManager.

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

            // Offset from the beginning of the graph. +1 because the beginning base will be at position 1.
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

        int getAltCount ();      // Return the number alternate alleles mapped to by the current alignment.
        std::string getCigarString (IAdjudicator::SharedPtr adjudicatorPtr);

        std::vector< uint32_t > getNodeIDs ();
        std::vector< int32_t > getNodeLengths ();
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
        // Overflow is occuring.
        int m_alt_count;     // Number of alternate alleles the alignment mapped to.
        
        // BED info.
        std::vector< uint32_t > m_node_IDs;
        std::vector< int32_t > m_node_lengths;
        std::vector< NodeInfo::VariantType > m_variant_types;

	};

}

#endif //GRAPHITE_GSSW_MAPPING_H
