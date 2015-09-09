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

	void GSSWAdjudicator::adjudicateMapping(IMapping::SharedPtr mappingPtr)
	{
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		auto swScore = mappingPtr->getMappingScore();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
		if (swPercent >= this->m_sw_percent)
		{
			auto incrementFunct = (alignmentPtr->isReverseStrand()) ? &IAllele::incrementReverseCount : &IAllele::incrementForwardCount;
			auto allelePtrs = mappingPtr->getAllelePtrs();
			for (uint32_t i = 0; i < allelePtrs.size(); ++i)
			{
				auto allelePtr = allelePtrs[i];
				auto mappingAlignmentInfoPtr = mappingPtr->getMappingAlignmentInfo(allelePtr, shared_from_this());
				auto alleleMappingScorePercent = ((mappingAlignmentInfoPtr->getSWScore() / (double)(mappingAlignmentInfoPtr->getLength() * this->m_match_value)) * 100);
				if (this->m_sw_percent > alleleMappingScorePercent) { continue; } // if the mapping score for this allele is too low then do not count it

				// check that the alignment maps to unique areas of the allele
				if (auto variantPtr = allelePtr->getVariantWPtr().lock()) // get the weakptr and lock it (checks if not expired)
				{
					if ((i == 0 && variantPtr->getAllelePrefixOverlapMaxCount(allelePtr) <= mappingAlignmentInfoPtr->getPrefixMatch()) || (i == allelePtrs.size() - 1 && variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr) <= mappingAlignmentInfoPtr->getSuffixMatch()))
					{
						continue;
					}
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
