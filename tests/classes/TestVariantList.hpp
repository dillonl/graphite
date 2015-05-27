#ifndef GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
#define GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP

#include "core/variant/VariantList.h"

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
		    TestVariantList(std::vector< Variant::SharedPtr >& variantPtrList)
			{
				m_current_index = 0;
				m_variants_ptr_list = std::vector< Variant::SharedPtr >(variantPtrList.begin(), variantPtrList.end());
			}

			~TestVariantList()
			{
			}

		};

	}
}

#endif //GWIZ_TESTS_CLASSES_TESTVARIANTLIST_HPP
