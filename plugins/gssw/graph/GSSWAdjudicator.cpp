#include "GSSWAdjudicator.h"
#include "GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "core/path/Path.h"
#include "core/mapping/MappingManager.h"
#include "AlignmentReporter.h"
#include "GSSWMapping.h"

#include <memory>

namespace gwiz
{
namespace gssw
{
	GSSWAdjudicator::GSSWAdjudicator(uint32_t swPercent, int matchValue) :
		m_sw_percent(swPercent),
		m_match_value(matchValue)
	{
	}

	GSSWAdjudicator::~GSSWAdjudicator()
	{
	}

	void GSSWAdjudicator::adjudicateMapping(IMapping::SharedPtr mappingPtr)
	{
		auto swScore = mappingPtr->getMappingScore();
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		float swPercent = (swScore / (alignmentPtr->getLength() * this->m_match_value));
		if (swPercent >= this->m_sw_percent)
		{
			if (alignmentPtr->isReverseStrand())
			{
				for (auto& allelePtr : mappingPtr->getAllelePtrs()) { allelePtr->incrementReverseCount(); }
			}
			else
			{
				for (auto& allelePtr : mappingPtr->getAllelePtrs()) { allelePtr->incrementForwardCount(); }
			}
		}
	}

}
}
