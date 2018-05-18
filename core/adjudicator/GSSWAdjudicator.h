#ifndef GRAPHITE_ADJUDICATOR_GSSWADJUDICATOR_H
#define GRAPHITE_ADJUDICATOR_GSSWADJUDICATOR_H

#include "core/util/Noncopyable.hpp"
#include "core/mapping/GSSWMapping.h"
#include "core/mapping/MappingAlignmentInfo.h"
#include "core/allele/IAllele.h"

#include <mutex>
#include <memory>

namespace graphite
{
	class GSSWMapping;
    class GSSWAdjudicator : private Noncopyable, public std::enable_shared_from_this< GSSWAdjudicator >
	{
	public:
		typedef std::shared_ptr< GSSWAdjudicator > SharedPtr;
		GSSWAdjudicator(uint32_t swPercent, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue);
		~GSSWAdjudicator();

		bool adjudicateMapping(std::shared_ptr< GSSWMapping > mappingPtr, uint32_t referenceSWPercent);
		int getMatchValue();
		int getMisMatchValue();
		int getGapOpenValue();
		int getGapExtensionValue();

		static uint32_t s_adj_count;

	private:
        void mapAllele(std::shared_ptr< IAllele > allelePtr, MappingAlignmentInfo::SharedPtr mappingAlignmentInfoPtr, std::shared_ptr< GSSWMapping > mappingPtr, IAlignment::SharedPtr alignmentPtr, bool referenceSWScoreIdentical, float swPercent);

		std::mutex m_adjudication_lock;
		uint32_t m_sw_percent;
		int m_match_value;
		int m_mismatch_value;
		int m_gap_open_value;
		int m_gap_extension_value;
	};
}

#endif //GRAPHITE_GSSW_GSSWADJUDICATOR_H
