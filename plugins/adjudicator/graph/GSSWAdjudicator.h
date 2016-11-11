#ifndef GRAPHITE_ADJUDICATOR_GSSWADJUDICATOR_H
#define GRAPHITE_ADJUDICATOR_GSSWADJUDICATOR_H

#include "core/adjudicator/IAdjudicator.h"

#include <mutex>

namespace graphite
{
namespace adjudicator
{
	class GSSWAdjudicator : public IAdjudicator
	{
	public:
		typedef std::shared_ptr< GSSWAdjudicator > SharedPtr;
		GSSWAdjudicator(uint32_t swPercent, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue);
		~GSSWAdjudicator();

		void adjudicateMapping(IMapping::SharedPtr mappingPtr) override;
		int getMatchValue() override;
		int getMisMatchValue() override;
		int getGapOpenValue() override;
		int getGapExtensionValue() override;

		static uint32_t s_adj_count;

	private:
		void mapAllele(IAllele::SharedPtr allelePtr, MappingAlignmentInfo::SharedPtr mappingAlignmentInfoPtr, IMapping::SharedPtr mappingPtr, IAlignment::SharedPtr alignmentPtr, bool checkAllelePrefix, bool checkAlleleSuffix);

		std::mutex m_adjudication_lock;
		uint32_t m_sw_percent;
		int m_match_value;
		int m_mismatch_value;
		int m_gap_open_value;
		int m_gap_extension_value;
	};
}
}

#endif //GRAPHITE_GSSW_GSSWADJUDICATOR_H
