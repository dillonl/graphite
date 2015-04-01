#ifndef GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
#define GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP

#include "core/variants/VariantList.h"

namespace gwiz
{
	namespace testing
	{

		/*
		 * This class is only for testing.
		 */

		class TestVariantList : public VariantList
		{
		public:
			typedef std::shared_ptr< TestVariantList > SharedPtr;
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
				bool hasVariants = (this->m_variant_index < this->m_variant_ptr_list.size());
				if (hasVariants)
				{
					variant = this->m_variant_ptr_list[this->m_variant_index];
					++this->m_variant_index;
				}
				return hasVariants;
			}

			void rewind() override
			{
				this->m_variant_index = 0;
			}

			size_t size()
			{
				return this->m_variant_ptr_list.size();
			}

			size_t getCount() override { return m_variant_ptr_list.size(); }

		private:

			size_t m_variant_index;
			std::vector< Variant::SharedPtr > m_variant_ptr_list;

		};

	}
}

#endif //GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
