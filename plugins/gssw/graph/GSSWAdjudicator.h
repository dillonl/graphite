#ifndef GWIZ_GSSW_GSSWADJUDICATOR_H
#define GWIZ_GSSW_GSSWADJUDICATOR_H

#include "core/adjudicator/IAdjudicator.h"

#include <mutex>

namespace gwiz
{
namespace gssw
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
	private:

		std::mutex m_adjudication_lock;
		uint32_t m_sw_percent;
		int m_match_value;
		int m_mismatch_value;
		int m_gap_open_value;
		int m_gap_extension_value;
	};
}
}

#endif //GWIZ_GSSW_GSSWADJUDICATOR_H
