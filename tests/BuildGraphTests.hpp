#include "gtest/gtest.h"

#include "TestConfig.h"

#include "core/graph/IGraph.h"
#include "core/reference/FastaReference.h"
#include "core/variants/Variant.h"
#include "core/variants/VariantList.h"
#include "core/variants/VCFParser.hpp"

namespace
{

// The fixture for testing class Foo.
	class BuildGraphTests : public ::testing::Test
	{
	protected:
		// You can remove any or all of the following functions if its body
		// is empty.

		BuildGraphTests()
		{
				// You can do set-up work for each test here.
		}

		virtual ~BuildGraphTests()
		{
			// You can do clean-up work that doesn't throw exceptions here.
	    }

		// If the constructor and destructor are not enough for setting up
		// and cleaning up each test, you can define the following methods:

		virtual void SetUp()
		{
			// Code here will be called immediately after the constructor (right
			// before each test).
		}

		virtual void TearDown()
		{
			// Code here will be called immediately after each test (right
			// before the destructor).
		}

		// Objects declared here can be used by all tests in the test case for Foo.
	};

	class TestGraphBGT : public gwiz::IGraph
	{
    public:
		TestGraphBGT(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantsListPtr) :
			gwiz::IGraph(referencePtr, variantsListPtr)
		{
		}
		~TestGraphBGT() {}
		void constructGraph() {}

		bool getNextCompoundVariantTest(gwiz::Variant::SharedPtr& variant)
		{
			return getNextCompoundVariant(variant);
		}

	};

	// Tests that the Foo::Bar() method does Abc.

	TEST_F(BuildGraphTests, CompoundVariantSimple)
	{
		auto regionPtr = std::make_shared< gwiz::Region >("20:50-80");
		std::string fastaPath = TEST_FASTA_FILE;
		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);

		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t60\t.\tTAAGT\tT\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t61\t.\tAA\tAGAAG\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t63\t.\tGTTAAC\tGTAG,G\t.\t.", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);
		variantsListPtr->addVariant(variantPtr3);

		auto testGraph = std::make_shared< TestGraphBGT >(fastaReferencePtr, variantsListPtr);
		gwiz::Variant::SharedPtr variantPtr;
		testGraph->getNextCompoundVariantTest(variantPtr);
		ASSERT_STREQ(variantPtr->getRef().c_str(), "TAAGTTAAC");
		ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "TTAAC");
		ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "TAGAAGGTTAAC");
		ASSERT_STREQ(variantPtr->getAlt()[2].c_str(), "TAAGTAG");
		ASSERT_STREQ(variantPtr->getAlt()[3].c_str(), "TAAG");
		ASSERT_EQ(variantPtr->getAlt().size(), 4);


		EXPECT_EQ(0, 0);
	}

	TEST_F(BuildGraphTests, BuildGraph)
	{
		auto regionPtr = std::make_shared< gwiz::Region >("20:20301000-20402000");
		std::string fastaPath = TEST_FASTA_FILE;
		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);

		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t20301046\t.\tTAATATATGTAATATATATTATATATGTAATATAATATATGTAAT\tT\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t20301055\t.\tTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT\tT\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t20301104\t.\tGT\tG\t.\t.", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);
		variantsListPtr->addVariant(variantPtr3);

		auto testGraph = std::make_shared< TestGraphBGT >(fastaReferencePtr, variantsListPtr);
		gwiz::Variant::SharedPtr variantPtr;
		testGraph->getNextCompoundVariantTest(variantPtr);

		ASSERT_STREQ(variantPtr->getRef().c_str(),    "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT");
		ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "TATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT");
		ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "TAATATATGT");
		ASSERT_STREQ(variantPtr->getAlt()[2].c_str(), "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT");
		ASSERT_EQ(variantPtr->getAlt().size(), 3);

		EXPECT_EQ(0, 0);
	}
}
