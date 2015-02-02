#ifndef GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
#define GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP

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
				if (this->m_variant_index == this->m_variant_ptr_list.size())
				{
					std::cout << "getnextvaraint [false]: " << this->m_variant_ptr_list.size()<< std::endl;
					return false;
				}
				else
				{
					std::cout << "getnextvaraint [true]: " << this->m_variant_ptr_list.size()<< std::endl;
					variant = this->m_variant_ptr_list[this->m_variant_index];
					++this->m_variant_index;
					return true;
				}
			}

			size_t size()
			{
				return this->m_variant_ptr_list.size();
			}

		private:
			size_t m_variant_index;
			std::vector< Variant::SharedPtr > m_variant_ptr_list;

		};

	}
}

#endif //GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
