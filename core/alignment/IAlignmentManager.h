#ifndef GRAPHITE_IALIGNMENTMANAGER_H
#define GRAPHITE_IALIGNMENTMANAGER_H

#include "AlignmentList.h"
#include "core/region/Region.h"
#include "core/util/Noncopyable.hpp"
#include "core/variant/IVariantList.h"
#include "core/variant/IVariantManager.h"
#include "core/sample/SampleManager.h"

#include <memory>
#include <functional>
#include <unordered_set>
#include <algorithm>

namespace graphite
{
	class Sample;
	class IAlignmentManager : private Noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignmentManager > SharedPtr;
	    IAlignmentManager() : m_name_alignment_ptr_map_ptr(std::make_shared< std::unordered_map< std::string, IAlignment::SharedPtr > >()) {}
		virtual ~IAlignmentManager() {}

		virtual void releaseResources() = 0;
		/* virtual void processMappingStatistics() = 0; */
		virtual SampleManager::SharedPtr getSamplePtrs() = 0;
		virtual IAlignmentList::SharedPtr getAlignmentsInRegion(Region::SharedPtr regionPtr)
		{
			position startPosition = regionPtr->getStartPosition();
			position endPosition = regionPtr->getEndPosition();

			auto lowerBound = std::lower_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [startPosition](const IAlignment::SharedPtr& alignmentPtr, const IAlignment::SharedPtr& ignore) {
					return startPosition > alignmentPtr->getPosition() + alignmentPtr->getLength();
				});
			auto upperBound = std::upper_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [endPosition](const IAlignment::SharedPtr& ignore, const IAlignment::SharedPtr& alignmentPtr) {
					return alignmentPtr->getPosition() > endPosition;
				});

			std::vector< IAlignment::SharedPtr > alignmentPtrs;
			alignmentPtrs.insert(alignmentPtrs.begin(), lowerBound, upperBound);
			return std::make_shared< AlignmentList >(alignmentPtrs,m_name_alignment_ptr_map_ptr, regionPtr);
		}

	protected:
		std::vector< Region::SharedPtr > getRegionsContainingVariantsWithPadding(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
		{
			std::vector< Region::SharedPtr > regionPtrs;

			auto variantListPtr = variantManagerPtr->getVariantsInRegion(this->m_region_ptr);
			std::vector< IVariant::SharedPtr > variantPtrs;
			IVariant::SharedPtr variantPtr;
			while (variantListPtr->getNextVariant(variantPtr))
			{
				variantPtrs.emplace_back(variantPtr);
			}

			Region::SharedPtr regionPtr = nullptr;
			for (auto i = 0; i < variantPtrs.size(); ++i)
			{
				uint32_t maxAlleleSize = 3000;
				position startPosition = ((variantPtrs[i]->getPosition() - (variantPadding * 2)) > 0) ? (variantPtrs[i]->getPosition() - (variantPadding * 2)) : 0;
				position endPosition = (variantPtrs[i]->getPosition() + maxAlleleSize + (variantPadding * 2));
				auto j = i + 1;
				while (j < (variantPtrs.size() - 1) && variantPtrs[j]->getPosition() < endPosition)
				{
					endPosition = (variantPtrs[j]->getPosition() + maxAlleleSize + (variantPadding * 2));
					++j;
				}
				i = j - 1;
				regionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), startPosition, endPosition, Region::BASED::ONE);
				regionPtrs.push_back(regionPtr);
			}
			return regionPtrs;
		}

		std::vector< IAlignment::SharedPtr > m_alignment_ptrs;
		std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > m_name_alignment_ptr_map_ptr;
		Region::SharedPtr m_region_ptr;
	};
}

#endif //GRAPHITE_IALIGNMENTMANAGER_H
