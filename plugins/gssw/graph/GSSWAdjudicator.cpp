#include "GSSWAdjudicator.h"
#include "GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "core/allele/EquivalentAllele.h"
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
		auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex(); // lock the alignmentMappingMutex so no one else can set the mapping ptr
		std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this
		if (auto alignmentMappingWPtr = alignmentPtr->getMapping().lock()) // check if the alignment has already been aligned previously
		{
			auto swAlignmentScore = alignmentMappingWPtr->getMappingScore();
			uint32_t swAlignmentPercent = ((swAlignmentScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
			if (swAlignmentPercent <= swPercent) { return; }
		}
		auto gsswMappingPtr = (GSSWMapping*)mappingPtr.get();
		auto mappingAlignmentInfoPtrs = gsswMappingPtr->getMappingAlignmentInfoPtrs(shared_from_this());
		// position remapPosition = gsswMappingPtr->position;
		bool mapped = false;
		std::unordered_map< IVariant::SharedPtr, bool > variantMap;
		auto incrementFunct = (alignmentPtr->isReverseStrand()) ? &IAllele::incrementReverseCount : &IAllele::incrementForwardCount;

		for (uint32_t i = 0; i < mappingAlignmentInfoPtrs.size(); ++i)
		{
			auto mappingAlignmentInfoPtr = mappingAlignmentInfoPtrs[i];
			auto allelePtr = mappingAlignmentInfoPtr->getAllelePtr();
			auto equivalentAllelePtr = std::dynamic_pointer_cast< EquivalentAllele >(allelePtr);
			std::vector< IAllele::SharedPtr > allelePtrs;
			if (equivalentAllelePtr)
			{
				allelePtrs = equivalentAllelePtr->getAllAlleles();
			}
			else
			{
				allelePtrs.emplace_back(allelePtr);
			}
			for (auto currentAllelePtr : allelePtrs)
			{
				auto variantPtr = currentAllelePtr->getVariantWPtr().lock();
				if (variantPtr == nullptr)
				{
					continue;
				}
				variantMap.emplace(variantPtr, true);
				if (swPercent < this->m_sw_percent) //if the percentage isn't high enough to increment the count
				{
					continue;
				}
				auto alleleMappingScorePercent = ((mappingAlignmentInfoPtr->getSWScore() / (double)(mappingAlignmentInfoPtr->getLength() * this->m_match_value)) * 100);
				if (this->m_sw_percent > alleleMappingScorePercent) // if the mapping score for this allele is too low then do not count it
				{
					continue;
				}

				// check that the alignment maps to unique areas of the allele
				if ((i == 0 && variantPtr->getAllelePrefixOverlapMaxCount(allelePtr) <= mappingAlignmentInfoPtr->getPrefixMatch()) || (i == mappingAlignmentInfoPtrs.size() - 1 && variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr) <= mappingAlignmentInfoPtr->getSuffixMatch()))
				{
					continue;
				}
				mapped = true;
				(*currentAllelePtr.*incrementFunct)();
			}
		}

		if (mapped)
		{
			alignmentPtr->setMapping(mappingPtr);
		}
		bool mappedToUnmapped = (!mapped && alignmentPtr->isMapped());
		bool unmappedToMapped = (mapped && !alignmentPtr->isMapped());

		if (mappedToUnmapped || unmappedToMapped)
		{
			for (auto iter : variantMap)
			{
				if (unmappedToMapped) { iter.first->incrementUnmappedToMappedCount(); }
				if (mappedToUnmapped) { iter.first->incrementMappedToUnmappedCount(); }
				if (mapped && alignmentPtr->getPosition() != gsswMappingPtr->getPosition()) { iter.first->incrementRepositionedCount(); }
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
