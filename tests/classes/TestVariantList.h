#ifndef GWIZ_TESTS_CLASSES_TESTVARIANTLIST_H
#define GWIZ_TESTS_CLASSES_TESTVARIANTLIST_H

#include "core/variants/IVariantList.h"

namespace gwiz
{
	namespace testing
	{

		/*
		 * This class is only for testing.
		 */

		class TestVariantList : public IVariantList
		{
		public:
		    TestVariantList(std::vector< Variant::SharedPtr > variantPtrList) :
			    m_variant_index(0),
			    m_variant_ptr_list(variantPtrList)
			{
			}

			~TestVariantList()
			{
			}

			bool getNextVariant(Variant::SharedPtr& variant) override
			{
				if (this->m_variant_index == this->m_variant_ptr_list.size())
				{
					return false;
				}
				else
				{
					variant = this->m_variant_ptr_list[this->m_variant_index];
					++this->m_variant_index;
					return true;
				}
			}

		private:
			size_t m_variant_index;
			std::vector< Variant::SharedPtr > m_variant_ptr_list;

		};

	}
}

#endif //GWIZ_TESTS_CLASSES_TESTVARIANTLIST_H
