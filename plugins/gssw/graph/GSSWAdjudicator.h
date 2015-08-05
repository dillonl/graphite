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
		GSSWAdjudicator(uint32_t swPercent, int matchValue);
		~GSSWAdjudicator();

		void adjudicateMapping(IMapping::SharedPtr mappingPtr) override;
	private:

		std::mutex m_adjudication_lock;
		uint32_t m_sw_percent;
		int m_match_value;
	};
}
}

#endif //GWIZ_GSSW_GSSWADJUDICATOR_H
