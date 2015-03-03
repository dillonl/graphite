#ifndef GWIZ_VARIANTLIST_H
#define GWIZ_VARIANTLIST_H

#include "IVariantList.h"

#include <memory>
#include <queue>

namespace gwiz
{
	class VariantList : public IVariantList
	{
	public:
		typedef std::shared_ptr< VariantList > SharedPtr;
		VariantList() {}
		~VariantList() {}

		inline bool getNextVariant(Variant::SharedPtr& variantPtr) override
		{
			bool hasVariants = !this->m_variants.empty();
			if (hasVariants)
			{
				variantPtr = this->m_variants.front();
				this->m_variants.pop();
			}
			return hasVariants;
		}

		inline void addVariant(const Variant::SharedPtr variant)
		{
			this->m_variants.push(variant);
		}

	private:
		std::queue< Variant::SharedPtr > m_variants;
	};
}

#endif //GWIZ_VARIANTLIST_H
