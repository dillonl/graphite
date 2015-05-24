#ifndef GWIZ_VARIANTLIST_H
#define GWIZ_VARIANTLIST_H

#include "IVariantList.h"

#include <typeinfo>
#include <memory>
#include <vector>
#include <algorithm>
#include <map>
#include <mutex>

namespace gwiz
{
	class VariantList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
	    VariantList() : m_current_index(0) {}
		virtual ~VariantList() {}

		inline bool getNextVariant(Variant::SharedPtr& variantPtr) override
		{
			std::lock_guard< std::mutex > lock(this->m_variants_ptr_list_mutex);
			bool hasVariants = this->m_current_index < this->m_variants_ptr_list.size();
			if (hasVariants)
			{
				variantPtr = this->m_variants_ptr_list[this->m_current_index];
				++this->m_current_index;
			}
			else
			{
				variantPtr = nullptr;
			}
			return hasVariants;
		}

		inline void addVariant(const Variant::SharedPtr variant) override
		{
			std::lock_guard< std::mutex > lock(this->m_variants_ptr_list_mutex);
			this->m_variants_ptr_list.emplace_back(variant);
		}

		inline void addVariants(IVariantList::SharedPtr variantsListPtr) override
		{
			auto variantsListPtrCast = std::dynamic_pointer_cast< VariantList >(variantsListPtr);
			if (variantsListPtrCast)
			{
				std::lock_guard< std::mutex > lock(this->m_variants_ptr_list_mutex);
				this->m_variants_ptr_list.insert(this->m_variants_ptr_list.end(), variantsListPtrCast->m_variants_ptr_list.begin(), variantsListPtrCast->m_variants_ptr_list.end());

			}
			else
			{
				Variant::SharedPtr variantPtr;
				while (variantsListPtr->getNextVariant(variantPtr))
				{
					addVariant(variantPtr);
				}
			}
		}

		void sort()
		{
			std::lock_guard< std::mutex > lock(this->m_variants_ptr_list_mutex);
			std::sort(this->m_variants_ptr_list.begin(), this->m_variants_ptr_list.end(), [](const Variant::SharedPtr& v1, const Variant::SharedPtr& v2) {
					return v1->getPosition() < v2->getPosition();
				});
		}

		inline void rewind()
		{
			m_current_index = 0;
		}

		IVariantList::SharedPtr getVariantsInRegion(Region::SharedPtr regionPtr) override
		{
			this->m_variants_ptr_list_mutex.lock();
			position startPosition = regionPtr->getStartPosition();
			position endPosition = regionPtr->getEndPosition();

			auto lowerBound = std::lower_bound(this->m_variants_ptr_list.begin(), this->m_variants_ptr_list.end(), nullptr, [startPosition](const Variant::SharedPtr& variantPtr, const Variant::SharedPtr& ignore) {
					return startPosition > variantPtr->getPosition();
				});
			auto upperBound = std::upper_bound(this->m_variants_ptr_list.begin(), this->m_variants_ptr_list.end(), nullptr, [endPosition](const Variant::SharedPtr& ignore, const Variant::SharedPtr& variantPtr) {
					return variantPtr->getPosition() > endPosition;
				});

			auto variantListPtr = std::make_shared< VariantList >();
			variantListPtr->m_variants_ptr_list = std::vector< Variant::SharedPtr >(lowerBound, upperBound);
			this->m_variants_ptr_list_mutex.unlock();
			return variantListPtr;
		}

		size_t getCount() override { return m_variants_ptr_list.size(); }

		void printVCFHeader(std::ostream& out)
		{
			out << "##fileformat=VCFv4.2" << std::endl;
			out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">" << std::endl;
			out << "##INFO=<ID=TC,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles, used in variant calling.\">" << std::endl;
			out << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth including non-counted low quality reads\">" << std::endl;
			out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
			out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" << std::endl;
		}

		void printToVCF(std::ostream& out) override
		{
			printVCFHeader(out);
			std::map< uint32_t, bool > variantPrintedMap;
			std::lock_guard< std::mutex > lock(this->m_variants_ptr_list_mutex);
			for (auto& variantPtr : this->m_variants_ptr_list)
			{
				if (variantPrintedMap.find(variantPtr->getVariantID()) != variantPrintedMap.end()) { continue; }
				/* if (!variantPtr->hasAlts()) { continue; } */
				variantPtr->printVariant(out);
				variantPrintedMap[variantPtr->getVariantID()] = true;
			}
		}

	protected:
		size_t m_current_index;
		std::vector< Variant::SharedPtr > m_variants_ptr_list;
		std::mutex m_variants_ptr_list_mutex;
	};
}

#endif //GWIZ_VARIANTLIST_H
