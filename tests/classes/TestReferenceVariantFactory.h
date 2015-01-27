#ifndef GWIZ_TEST_TESTREFERENCEVARIANTFACTORY_H
#define GWIZ_TEST_TESTREFERENCEVARIANTFACTORY_H

namespace gwiz
{
	namespace test
	{
		class TestReferenceVariantFactory
		{
		public:
			static void CreateVariants()
			{
			}

			static void (std::string& reference, std::vector< std::string > vcf_list,
						 Region::SharedPtr region, IReference::SharedPtr reference,
						 IVariantList::SharedPtr variantList)
			{

			}

		private:
			TestReferenceVariantFactory() {}
			~TestReferenceVariantFactory() {}
		};
	}

}

#endif //GWIZ_TEST_TESTREFERENCEVARIANTFACTORY_H
