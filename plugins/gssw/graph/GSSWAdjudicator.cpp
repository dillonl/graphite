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
		// static std::mutex s_lock;
		// std::lock_guard< std::mutex > sGuard(s_lock);


		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex(); // lock the alignmentMappingMutex so no one else can set the mapping ptr
		std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this

		auto swScore = mappingPtr->getMappingScore();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);

		// std::cout << "mapping id: " <<  mappingPtr->m_id << std::endl;
		// std::cout << "swPercent: " << swPercent << std::endl;
		if (swPercent < this->m_sw_percent) //if the percentage isn't high enough to increment the count
		{
			// std::cout << "killed on percent: " << swPercent << " " << this->m_sw_percent << std::endl;
			return;
		}
		auto alignmentMappingWPtr = alignmentPtr->getMapping().lock();
		if (alignmentMappingWPtr && alignmentMappingWPtr->getMapped()) // check if the alignment has already been aligned previously
		{
			auto swAlignmentScore = alignmentMappingWPtr->getMappingScore();
			uint32_t swAlignmentPercent = ((swAlignmentScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
			if (swAlignmentPercent <= swPercent) // if the new alignment isn't "better" than the previous alignment
			{
				// std::cout << "killed on previous mapping percent: " << swAlignmentPercent << std::endl;
				return;
			}
			alignmentMappingWPtr->setMapped(false);
		}

		auto gsswMappingPtr = (GSSWMapping*)mappingPtr.get();
		auto mappingAlignmentInfoPtrs = gsswMappingPtr->getMappingAlignmentInfoPtrs(shared_from_this());
		bool mapped = false;

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
				auto alleleMappingScorePercent = ((mappingAlignmentInfoPtr->getSWScore() / (double)(mappingAlignmentInfoPtr->getLength() * this->m_match_value)) * 100);
				if ((variantPtr == nullptr) || // if this allele doesn't belong to a variant in the VCF (if it's a reference node at a non variant site)
					(alleleMappingScorePercent < this->m_sw_percent) || // if the mapping score for this allele is too low then do not count it
					((i == 0 && variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr) >= mappingAlignmentInfoPtr->getPrefixMatch()) || // check that the alignment maps to unique areas of the allele
					 (i == mappingAlignmentInfoPtrs.size() - 1 && variantPtr->getAllelePrefixOverlapMaxCount(allelePtr) >= mappingAlignmentInfoPtr->getSuffixMatch())))
				{
					continue;
				}
				else // if this allele was successfully mapped,
				{
					mapped = true;
					auto countFunct = (alignmentPtr->isReverseStrand())? &IAllele::incrementReverseCount : &IAllele::incrementForwardCount;
					mappingPtr->addAlleleCountCallback(std::bind(countFunct, currentAllelePtr));
					alignmentPtr->setMapping(mappingPtr);
					mappingPtr->setMapped(true);
				}
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
