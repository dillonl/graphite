#include "GSSWAdjudicator.h"
#include "core/graph/GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "core/allele/EquivalentAllele.h"
#include "core/mapping/MappingManager.h"
#include "core/mapping/GSSWMapping.h"

#include <memory>

namespace graphite
{
	uint32_t GSSWAdjudicator::s_adj_count = 0;
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

	bool GSSWAdjudicator::adjudicateMapping(IMapping::SharedPtr mappingPtr, uint32_t referenceSWPercent)
	{
		// std::cout << "comment out this lock" << std::endl;
		// std::lock_guard< std::mutex > sGuard(s_lock);
		// {
			// std::lock_guard< std::mutex > sGuard(s_lock);
			// ++s_adj_count;
		// }
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		if (alignmentPtr->getLength() <= 0 || alignmentPtr->getSequence() == nullptr)
		{
			std::cout << "no sequence: " << alignmentPtr->getPosition() << std::endl;
		}
		auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex(); // lock the alignmentMappingMutex so no one else can set the mapping ptr
		std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this

		auto swScore = mappingPtr->getMappingScore();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);

		auto alignmentMappingWPtr = alignmentPtr->getMapping().lock();
		if (alignmentMappingWPtr && alignmentMappingWPtr->getMapped()) // check if the alignment has already been aligned previously
		{
			auto swAlignmentScore = alignmentMappingWPtr->getMappingScore();
			uint32_t swAlignmentPercent = ((swAlignmentScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);
			if (swAlignmentPercent <= swPercent) // if the new alignment isn't "better" than the previous alignment
			{
				return false;
			}
			alignmentMappingWPtr->setMapped(false);
		}

		auto gsswMappingPtr = (GSSWMapping*)mappingPtr.get();
		auto mappingAlignmentInfoPtrs = gsswMappingPtr->getMappingAlignmentInfoPtrs(shared_from_this());
		bool referenceSWScoreIdentical = swPercent == referenceSWPercent;

		// static std::mutex tmpLock;
		// std::lock_guard< std::mutex > l(tmpLock);
	    // std::cout << alignmentPtr->getPosition() << " " << alignmentPtr->getID() << " sw: " << swPercent << " refSW: " << referenceSWPercent << std::endl;

		bool didMap = false;

		for (uint32_t i = 0; i < mappingAlignmentInfoPtrs.size(); ++i)
		{
			auto mappingAlignmentInfoPtr = mappingAlignmentInfoPtrs[i];
			auto allelePtr = mappingAlignmentInfoPtr->getAllelePtr();
			auto equivalentAllelePtr = std::dynamic_pointer_cast< EquivalentAllele >(allelePtr);
			bool checkAllelePrefix = (i == mappingAlignmentInfoPtrs.size() - 1);
			bool checkAlleleSuffix = (i == 0);
			if (equivalentAllelePtr)
			{
				for (auto currentAllelePtr : equivalentAllelePtr->getAllAlleles())
				{
					didMap = didMap || mapAllele(currentAllelePtr, mappingAlignmentInfoPtr, mappingPtr, alignmentPtr, checkAllelePrefix, checkAlleleSuffix, referenceSWScoreIdentical);
				}
			}
			else
			{
				didMap = didMap || mapAllele(allelePtr, mappingAlignmentInfoPtr, mappingPtr, alignmentPtr, checkAllelePrefix, checkAlleleSuffix, referenceSWScoreIdentical);
			}

		}
		return didMap;
	}

	bool GSSWAdjudicator::mapAllele(IAllele::SharedPtr allelePtr, MappingAlignmentInfo::SharedPtr mappingAlignmentInfoPtr, IMapping::SharedPtr mappingPtr, IAlignment::SharedPtr alignmentPtr, bool checkAllelePrefix, bool checkAlleleSuffix, bool referenceSWScoreIdentical)
	{
		auto variantPtr = allelePtr->getVariantWPtr().lock();
		auto alleleMappingScorePercent = ((mappingAlignmentInfoPtr->getSWScore() / (double)(mappingAlignmentInfoPtr->getLength() * this->m_match_value)) * 100);

		if (variantPtr != nullptr) // || // if this allele doesn't belong to a variant in the VCF (if it's a reference node at a non variant site)
		{
			bool isAlt = allelePtr->getID() % 2 != 0;
			// if ((referenceSWScoreIdentical && isAlt) || (checkAlleleSuffix && variantPtr->getAlleleSuffixOverlapMaxCount(allelePtr) >= mappingAlignmentInfoPtr->getPrefixMatch()) ||(checkAllelePrefix && variantPtr->getAllelePrefixOverlapMaxCount(allelePtr) >= mappingAlignmentInfoPtr->getSuffixMatch()))  // check that the alignment maps to unique areas of the allele
			if (referenceSWScoreIdentical && isAlt)
			{
				mappingPtr->addAlleleCountCallback(std::bind(&IAllele::incrementCount, allelePtr, alignmentPtr, AlleleCountType::Ambiguous));
				alignmentPtr->setMapping(mappingPtr);
				mappingPtr->setMapped(true);
				/*
				{
					static std::mutex l;
					std::lock_guard< std::mutex > lock(l);
					mappingPtr->printMapping();
					// std::cout << mappingAlignmentInfoPtr->getSWScore() << std::endl;
				}
				*/
				return true;
			}
			else
			{
				mappingPtr->addAlleleCountCallback(std::bind(&IAllele::incrementCount, allelePtr, alignmentPtr, scoreToAlleleCountType(alleleMappingScorePercent)));
				alignmentPtr->setMapping(mappingPtr);
				mappingPtr->setMapped(true);
				return true;
			}
		}
		return false;
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
