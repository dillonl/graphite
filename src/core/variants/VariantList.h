#ifndef GWIZ_VARIANTLIST_H
#define GWIZ_VARIANTLIST_H

#include "IVariantList.h"

#include <memory>
#include <vector>

namespace gwiz
{
	class VariantList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
	    VariantList() : m_current_index(0) {}
		~VariantList() {}

		inline bool getNextVariant(Variant::SharedPtr& variantPtr) override
		{
			bool hasVariants = this->m_current_index < this->m_variants.size();
			if (hasVariants)
			{
				variantPtr = this->m_variants[this->m_current_index];
				++this->m_current_index;
			}
			return hasVariants;
		}

		inline void addVariant(const Variant::SharedPtr variant)
		{
			this->m_variants.emplace_back(variant);
		}

		inline void rewind()
		{
			m_current_index = 0;
		}

	private:
		size_t m_current_index;
		std::vector< Variant::SharedPtr > m_variants;
	};
}

#endif //GWIZ_VARIANTLIST_H
