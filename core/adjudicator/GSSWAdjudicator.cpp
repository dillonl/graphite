#include "GSSWAdjudicator.h"
#include "core/graph/GSSWGraph.h"
#include "core/variant/VariantList.h"
#include "core/allele/EquivalentAllele.h"
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

	bool GSSWAdjudicator::adjudicateMapping(std::shared_ptr< GSSWMapping > mappingPtr, uint32_t referenceSWPercent)
	{
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		if (alignmentPtr->getLength() <= 0 || alignmentPtr->getSequence() == nullptr)
		{
			std::cout << "no sequence: " << alignmentPtr->getPosition() << std::endl;
		}
		auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex(); // lock the alignmentMappingMutex so no one else can set the mapping ptr
		std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this

		auto swScore = mappingPtr->getMappingScore();
		uint32_t swPercent = ((swScore / (double)(alignmentPtr->getLength() * this->m_match_value)) * 100);

		auto gsswMappingPtr = (GSSWMapping*)mappingPtr.get();
		auto mappingAlignmentInfoPtrs = gsswMappingPtr->getMappingAlignmentInfoPtrs(shared_from_this());
		bool referenceSWScoreIdentical = swPercent == referenceSWPercent;
		bool didMap = false;

		for (uint32_t i = 0; i < mappingAlignmentInfoPtrs.size(); ++i)
		{
			auto mappingAlignmentInfoPtr = mappingAlignmentInfoPtrs[i];
			auto allelePtr = mappingAlignmentInfoPtr->getAllelePtr();
			auto equivalentAllelePtr = std::dynamic_pointer_cast< EquivalentAllele >(allelePtr);
			std::vector< IAllele::SharedPtr > allelePtrs = (equivalentAllelePtr) ? equivalentAllelePtr->getAllAlleles() : (std::vector< IAllele::SharedPtr >){allelePtr};
			for (auto currentAllelePtr : allelePtrs)
			{
				auto variantPtr = allelePtr->getVariantWPtr().lock();
				auto alleleMappingScorePercent = ((mappingAlignmentInfoPtr->getSWScore() / (double)(mappingAlignmentInfoPtr->getLength() * this->m_match_value)) * 100);
				if (variantPtr != nullptr && alleleMappingScorePercent > 0) // || // if this allele doesn't belong to a variant in the VCF (if it's a reference node at a non variant site) and if this matches in any way
				{
					mapAllele(currentAllelePtr, mappingAlignmentInfoPtr, mappingPtr, alignmentPtr, referenceSWScoreIdentical, swPercent);
					didMap = true;
				}
			}

		}
		return didMap;
	}

	void GSSWAdjudicator::mapAllele(IAllele::SharedPtr allelePtr, MappingAlignmentInfo::SharedPtr mappingAlignmentInfoPtr, GSSWMapping::SharedPtr mappingPtr, IAlignment::SharedPtr alignmentPtr, bool referenceSWScoreIdentical, float swPercent)
	{
		static int count = 0;
		bool isAlt = allelePtr->getID() % 2 != 0;
		AlleleCountType alleleCountType = (referenceSWScoreIdentical && isAlt) ? AlleleCountType::Ambiguous : scoreToAlleleCountType(swPercent);
		allelePtr->incrementCount(alignmentPtr->isReverseStrand(), alignmentPtr->getSample(), alleleCountType);
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
