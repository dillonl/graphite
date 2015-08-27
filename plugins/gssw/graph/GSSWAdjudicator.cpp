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
				auto mappingAlignInfoPtr = mappingPtr->getGSSWAlignmentPtrFromAllelePtr(allelePtr);
				std::cout << "Mapping: " << mappingAlignInfoPtr->getMappingOffset() << " " << mappingAlignInfoPtr->getMappingSequence().size() << " " <<  alignmentPtr->getLength() << std::endl;
				if (mappingAlignInfoPtr->getMappingOffset() == 0 && memcmp(allelePtr->getSequence(), mappingAlignInfoPtr->getMappingSequence().c_str(), mappingAlignInfoPtr->getMappingSequence().size()) == 0) // if the mapping is at the beginning then check the prefix for differences
				{
					continue;
				}
				if (memcmp((allelePtr->getSequence() + allelePtr->getLength() - mappingAlignInfoPtr->getMappingSequence().size()), mappingAlignInfoPtr->getMappingSequence().c_str(), mappingAlignInfoPtr->getMappingSequence().size()) == 0)
				{
					std::cout << "Mapped to the END!!!" << std::endl;
					continue;
				}
				// else if ((mappingAlignInfoPtr->getMappingOffset() + ) == alignmentPtr->getLength()) // if the mapping is at the end check the suffix
				// {
				// std::cout << "Mapped to the END!!!" << std::endl;
					// continue;
				// }
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
