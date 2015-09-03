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
		static std::mutex mlock;
		std::lock_guard< std::mutex > guard(mlock);
		std::cout << "locked" << std::endl;

		// ((GSSWMapping*)mappingPtr.get())->printLongFormat();

		auto swScore = mappingPtr->getMappingScore();
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
		if (swPercent >= this->m_sw_percent)
		{
			auto incrementFunct = (alignmentPtr->isReverseStrand()) ? &IAllele::incrementReverseCount : &IAllele::incrementForwardCount;
			uint32_t count = 0;
			for (auto& allelePtr : mappingPtr->getAllelePtrs())
			{
				uint32_t alleleMappingScorePercent = ((mappingPtr->getMappingScore() / (double)(allelePtr->getLength() * this->m_match_value)) * 100);
				/*
				std::cout << "ms: " << mappingPtr->getMappingScore() << std::endl;
				std::cout << "al: " << allelePtr->getLength() << std::endl;
				std::cout << "mv: " << m_match_value << std::endl;
				std::cout << "allele score: " << alleleMappingScorePercent << std::endl;
				*/
				if (this->m_sw_percent > alleleMappingScorePercent) { continue; } // if the mapping score for this allele is too low then do not count it

				// check that the alignment maps to unique areas of the allele
				if (auto variantPtr = allelePtr->getVariantWPtr().lock()) // get the weakptr and lock it (checks if not expired)
				{
					auto mappingAlignInfoPtr = mappingPtr->getGSSWAlignmentPtrFromAllelePtr(allelePtr);
					std::cout << "mapping length: " << mappingAlignInfoPtr->getMappingOffsetLength() << std::endl;
					// check the mapping prefix and suffix to make sure the alignment covers more than the pre/suffix
					if ((count == 0 && mappingAlignInfoPtr->getMappingOffsetLength() <= variantPtr->getAllelePrefixOverlapMaxCount(allelePtr)) ||
						(count == (mappingPtr->getAllelePtrs().size() - 1) && mappingAlignInfoPtr->getMappingOffsetLength() <= variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr)))
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
