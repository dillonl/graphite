#include "GSSWAdjudicator.h"
#include "GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "core/path/Path.h"
#include "core/mapping/MappingManager.h"
#include "AlignmentReporter.h"
#include "GSSWMapping.h"

#include <memory>

namespace graphite
{
namespace gssw
{
	GSSWAdjudicator::GSSWAdjudicator(uint32_t swPercent, int matchValue, int misMatchValue, int gapOpenValue, int gapExtensionValue) :
		m_sw_percent(swPercent),
		m_match_value(matchValue),
		m_mismatch_value(misMatchValue),
		m_gap_open_value(gapOpenValue),
		m_gap_extension_value(gapExtensionValue)
	{
	}

	GSSWAdjudicator::~GSSWAdjudicator()
	{
	}

	/*
	  typedef struct {
	  uint16_t score1;
	  uint16_t score2;
	  int32_t ref_begin1;
	  int32_t ref_end1;
	  int32_t	read_begin1;
	  int32_t read_end1;
	  int32_t ref_end2;
	  gssw_seed seed;
	  uint8_t is_byte;
	  void* mH;
	  } gssw_align;
	*/

	void GSSWAdjudicator::adjudicateMapping(IMapping::SharedPtr mappingPtr)
	{
		auto swScore = mappingPtr->getMappingScore();
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
		if (swPercent >= this->m_sw_percent)
		{
			auto incrementFunct = (alignmentPtr->isReverseStrand()) ? &IAllele::incrementReverseCount : &IAllele::incrementForwardCount;
			for (auto& allelePtr : mappingPtr->getAllelePtrs())
			{
				if (auto variantPtr = allelePtr->getVariantWPtr().lock()) // get the weakptr and lock it (checks if not expired)
				{
					auto prefixOverlapCount = variantPtr->getAllelePrefixOverlapMaxCount(allelePtr);
					auto suffixOverlapCount = variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr);

					auto mappingAlignInfoPtr = mappingPtr->getGSSWAlignmentPtrFromAllelePtr(allelePtr);
					/*
					if (mappingOffsetLength <= prefixOverlapCount && memcmp(allelePtr->getSequence(), mappingAlignInfoPtr->getMappingSequence().c_str(), prefixOverlapCount) == 0)
					{
						continue;
					}
					if (mappingOffsetLength <= suffixOverlapCount && memcmp(allelePtr->getSequence(), mappingAlignInfoPtr->getMappingSequence().c_str(), prefixOverlapCount) == 0)
					{
					}
					*/
				}

				(*allelePtr.*incrementFunct)();
			}
		}
	}

	int GSSWAdjudicator::getMatchValue()
	{
		return this->m_match_value;
	}

	int GSSWAdjudicator::getMisMatchValue()
	{
		return this->m_mismatch_value;
	}

	int GSSWAdjudicator::getGapOpenValue()
	{
		return this->m_gap_open_value;
	}

	int GSSWAdjudicator::getGapExtensionValue()
	{
		return this->m_gap_extension_value;
	}

}
}
