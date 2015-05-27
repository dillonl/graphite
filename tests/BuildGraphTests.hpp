#include "gtest/gtest.h"

#include <algorithm>

#include "TestConfig.h"

#include "core/graph/IGraph.h"
#include "core/reference/FastaReference.h"
#include "core/variant/Variant.h"
#include "core/variant/VariantList.h"
#include "core/variant/VCFParser.hpp"

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

	};

	// Tests that the Foo::Bar() method does Abc.

	TEST_F(BuildGraphTests, CompoundVariantSimple)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t60\t.\tTAAGT\tT\t.\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t61\t.\tAA\tAGAAG\t.\t.\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t63\t.\tGTTAAC\tGTAG,G\t.\t.\t.", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);
		variantsListPtr->addVariant(variantPtr3);

		gwiz::Variant::SharedPtr variantPtr;
		variantsListPtr->getNextCompoundVariant(variantPtr);
		ASSERT_STREQ(variantPtr->getRef().c_str(), "TAAGTTAAC");
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TTAAC") != variantPtr->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TAGAAGGTTAAC") != variantPtr->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TAAGTAG") != variantPtr->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TAAG") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtr->getAlt().size(), 4);
	}

	TEST_F(BuildGraphTests, BuildGraph)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t20301046\t.\tTAATATATGTAATATATATTATATATGTAATATAATATATGTAAT\tT\t.\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t20301055\t.\tTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT\tT\t.\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t20301104\t.\tGT\tG\t.\t.\t.", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);
		variantsListPtr->addVariant(variantPtr3);

		gwiz::Variant::SharedPtr variantPtr;
		variantsListPtr->getNextCompoundVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRef().c_str(),    "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT");
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT") != variantPtr->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TAATATATGT") != variantPtr->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtr->getAlt().size(), 3);
	}

	TEST_F(BuildGraphTests, BuildGraph2)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("1\t105774\t.\tGT\tG\t108.73\ttAC=1;AF=0.500;AN=2;BaseQRankSum=-1.532;ClippingRankSum=0.143;DP=11;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.91;MQ0=0;MQRankSum=0.466;QD=9.88;ReadPosRankSum=0.143\tGT:AD:DP:GQ:PL\t0/1:7,5:12:99:167,0,249\t", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("1\t105774\t.\tGTT\tGT\t3.22758\t.\tAB=0.277778;ABP=10.7311;AC=1;AF=0.5;AN=2;AO=5;CIGAR=1M1D1M;DP=18;DPB=16.3333;DPRA=0;EPP=6.91895;EPPR=3.17734;GTI=0;LEN=1;MEANALT=1;MQM=21.8;MQMR=16.5385;NS=1;NUMALT=1;ODDS=0.0976763;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=185;QR=465;RO=13;RPP=6.91895;RPPR=3.17734;RUN=1;SAF=5;SAP=13.8677;SAR=0;SRF=13", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);

		gwiz::Variant::SharedPtr variantPtr;
		variantsListPtr->getNextCompoundVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRef().c_str(),    "GTT");
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "GT") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtr->getAlt().size(), 1);
	}

	TEST_F(BuildGraphTests, BuildGraph3)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t69506\t.\tG\tGACAC\t590.52\t.\tAC=2;AF=1.00;AN=2;DP=41;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=3.60\tGT:AD:DP:GQ:PL\t1/1:0,25:25:82:1163,82,0", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t69506\t.\tGACACACACACACACACACACACACACACA\tGACACACACACACACACACACACACACACACACA\t506.892\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=22;CIGAR=1M4I29M;DP=30;DPB=46.5333;DP", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);

		gwiz::Variant::SharedPtr variantPtr;
		gwiz::Variant::SharedPtr variantPtrTest;
		variantsListPtr->getNextCompoundVariant(variantPtr);
		variantsListPtr->getNextCompoundVariant(variantPtrTest);

		ASSERT_STREQ(variantPtr->getRef().c_str(),    "GACACACACACACACACACACACACACACA");
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "GACACACACACACACACACACACACACACACACA") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtr->getAlt().size(), 1);
		ASSERT_EQ(variantPtrTest, nullptr);
	}

	TEST_F(BuildGraphTests, BuildGraph4)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t72104\t.\tTA\tT\t1205.73\t.\tAC=2;AF=1.00;AN=2;DP=39;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.92\tGT:AD:DP:GQ:PL  1/1:0,36:36:99:1243,109,0", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t72104\t.\tTAA\tTA\t1062.05\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=36;CIGAR=1M1D1M;DP=39;DPB=59.6667;DPRA=0;EPP=9.04217;EPPR=5.18177;GTI=0;LEN=1;MEANALT=3;MQM=6", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t72719\t.\tC\tT\t1779.36\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.63072;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantPtr1);
		variantsListPtr->addVariant(variantPtr2);
		variantsListPtr->addVariant(variantPtr3);

		gwiz::Variant::SharedPtr variantPtr;
		gwiz::Variant::SharedPtr variantPtrTest;
		gwiz::Variant::SharedPtr variantPtrNULL;
		variantsListPtr->getNextCompoundVariant(variantPtr);
		variantsListPtr->getNextCompoundVariant(variantPtrTest);
		variantsListPtr->getNextCompoundVariant(variantPtrNULL);

		ASSERT_STREQ(variantPtr->getRef().c_str(),    "TAA");
		ASSERT_TRUE(std::find(variantPtr->getAlt().begin(), variantPtr->getAlt().end(), "TA") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtr->getAlt().size(), 1);
		ASSERT_EQ(variantPtr->getPosition(), 72104);

		ASSERT_STREQ(variantPtrTest->getRef().c_str(),    "C");
		ASSERT_TRUE(std::find(variantPtrTest->getAlt().begin(), variantPtr->getAlt().end(), "T") != variantPtr->getAlt().end());
		ASSERT_EQ(variantPtrTest->getAlt().size(), 1);
		ASSERT_EQ(variantPtrTest->getPosition(), 72719);

		ASSERT_EQ(variantPtrNULL, nullptr);
	}

	TEST_F(BuildGraphTests, BuildGraph5)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantLinePtr1 = gwiz::Variant::BuildVariant("20\t86005\t.\tC\tCA\t31.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=-0.361;ClippingRankSum=0.361;DP=41;FS=7.782;MLEAC=1;MLEAF=0.500;MQ=58.77;MQ0=0;MQRankSum=0", vcfParser);
		auto variantLinePtr2 = gwiz::Variant::BuildVariant("20\t86164\t.\tAG\tA\t2371.73\t.\tAC=2;AF=1.00;AN=2;DP=65;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.22     GT:AD:DP:GQ:PL  1/1:0,65:65:99:2409,196,0", vcfParser);
		auto variantLinePtr3 = gwiz::Variant::BuildVariant("20\t86164\t.\tAGG\tAG\t2037.92\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=65;CIGAR=1M1D1M;DP=66;DPB=45;DPRA=0;EPP=3.04371;EPPR=0;GTI=0;LEN=1;MEANALT=2;MQM=59.8615;MQMR", vcfParser);
		auto variantLinePtr4 = gwiz::Variant::BuildVariant("20\t87263\t.\tGA\tG\t1556.73\t.\tAC=2;AF=1.00;AN=2;DP=60;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=25.95     GT:AD:DP:GQ:PL  1/1:0,53:53:99:1594,160,0", vcfParser);
		auto variantLinePtr5 = gwiz::Variant::BuildVariant("20\t87263\t.\tGAAAAAAAAAT\tGAAAAAAAAT\t1679.52\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=54;CIGAR=1M1D9M;DP=55;DPB=56.4545;DPRA=0;EPP=3.6537;EPPR=0;GTI=0;LEN=1;MEANAL", vcfParser);
		auto variantLinePtr6 = gwiz::Variant::BuildVariant("20\t87416\t.\tA\tC\t1752.74\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.16541;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N", vcfParser);

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantLinePtr1);
		variantsListPtr->addVariant(variantLinePtr2);
		variantsListPtr->addVariant(variantLinePtr3);
		variantsListPtr->addVariant(variantLinePtr4);
		variantsListPtr->addVariant(variantLinePtr5);
		variantsListPtr->addVariant(variantLinePtr6);

		gwiz::Variant::SharedPtr variantPtr1;
		gwiz::Variant::SharedPtr variantPtr2;
		gwiz::Variant::SharedPtr variantPtr3;
		gwiz::Variant::SharedPtr variantPtr4;
		gwiz::Variant::SharedPtr variantPtr5;
		variantsListPtr->getNextCompoundVariant(variantPtr1);
		variantsListPtr->getNextCompoundVariant(variantPtr2);
		variantsListPtr->getNextCompoundVariant(variantPtr3);
		variantsListPtr->getNextCompoundVariant(variantPtr4);
		variantsListPtr->getNextCompoundVariant(variantPtr5);

		ASSERT_EQ(variantPtr1->getPosition(), 86005);
		ASSERT_STREQ(variantPtr1->getRef().c_str(),    "C");
		ASSERT_TRUE(std::find(variantPtr1->getAlt().begin(), variantPtr1->getAlt().end(), "CA") != variantPtr1->getAlt().end());
		ASSERT_EQ(variantPtr1->getAlt().size(), 1);

		ASSERT_EQ(variantPtr2->getPosition(), 86164);
		ASSERT_STREQ(variantPtr2->getRef().c_str(),    "AGG");
		ASSERT_TRUE(std::find(variantPtr2->getAlt().begin(), variantPtr2->getAlt().end(), "AG") != variantPtr2->getAlt().end());
		ASSERT_EQ(variantPtr2->getAlt().size(), 1);

		ASSERT_EQ(variantPtr3->getPosition(), 87263);
		ASSERT_STREQ(variantPtr3->getRef().c_str(),    "GAAAAAAAAAT");
		ASSERT_TRUE(std::find(variantPtr3->getAlt().begin(), variantPtr3->getAlt().end(), "GAAAAAAAAT") != variantPtr3->getAlt().end());
		ASSERT_EQ(variantPtr3->getAlt().size(), 1);

		ASSERT_EQ(variantPtr4->getPosition(), 87416);
		ASSERT_STREQ(variantPtr4->getRef().c_str(),    "A");
		ASSERT_TRUE(std::find(variantPtr4->getAlt().begin(), variantPtr4->getAlt().end(), "C") != variantPtr4->getAlt().end());
		ASSERT_EQ(variantPtr4->getAlt().size(), 1);

		ASSERT_EQ(variantPtr5, nullptr);
	}

	TEST_F(BuildGraphTests, BuildGraph6)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantLinePtr1 = gwiz::Variant::BuildVariant("1\t249240051\t.\tTG\tTAG\t773.046\t.\tAB=0.506944;ABP=3.07062;AC=1;AF=0.5;AN=2;AO=73;CIGAR=1M1I1M;DP=144;DPB=191;DPRA=0;EPP=6.60959;EPPR=3.21711;GTI=0;LEN=1;MEANALT=8;MQM=21.7534", vcfParser);
		auto variantLinePtr2 = gwiz::Variant::BuildVariant("2\t249240099\t.\tTG\tTAG,TAGG\t151.835\t.\tAB=0.238636,0.227273;ABP=55.2243,59.8634;AC=1,1;AF=0.5,0.5;AN=2;AO=21,20;CIGAR=1M1I1M,1M2I1M;DP=88;DPB=149.5;DPRA=0,0;EPP=3.1137,3.0", vcfParser);
		/*
		auto variantLinePtr3 = gwiz::Variant::BuildVariant("1\t249240303\t.\tT\tTA\t1045.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=0.552;ClippingRankSum=-1.541;DP=42;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=34.41;MQ0=0;MQRankSum=-1.465;QD=24.90;Rea", vcfParser);
		auto variantLinePtr4 = gwiz::Variant::BuildVariant("1\t249240303\t.\tTG\tTAG\t654.301\t.\tAB=0.309028;ABP=94.2423;AC=1;AF=0.5;AN=2;AO=89;CIGAR=1M1I1M;DP=288;DPB=403.5;DPRA=0;EPP=8.49998;EPPR=8.28388;GTI=0;LEN=1;MEANALT=13;MQM=19.4", vcfParser);
		auto variantLinePtr5 = gwiz::Variant::BuildVariant("1\t249240537\t.\tGGTG\tGG,GGGG\t66.6045\t.\tAB=0.382353,0.264706;ABP=7.09778,19.3602;AC=1,1;AF=0.5,0.5;AN=2;AO=13,9;CIGAR=1M2D1M,2M1X1M;DP=34;DPB=29.75;DPRA=0,0;EPP=11.1951,3.25157;EPP", vcfParser);
		*/

		auto variantsListPtr = std::make_shared< gwiz::VariantList >();
		variantsListPtr->addVariant(variantLinePtr1);
		variantsListPtr->addVariant(variantLinePtr2);
		/*
		variantsListPtr->addVariant(variantLinePtr3);
		variantsListPtr->addVariant(variantLinePtr4);
		variantsListPtr->addVariant(variantLinePtr5);
		*/

		gwiz::Variant::SharedPtr variantPtr1;
		gwiz::Variant::SharedPtr variantPtr2;
		/*
		gwiz::Variant::SharedPtr variantPtr3;
		gwiz::Variant::SharedPtr variantPtr4;
		gwiz::Variant::SharedPtr variantPtr5;
		*/
		variantsListPtr->getNextCompoundVariant(variantPtr1);
		variantsListPtr->getNextCompoundVariant(variantPtr2);
		/*
		variantsListPtr->getNextCompoundVariant(variantPtr3);
		variantsListPtr->getNextCompoundVariant(variantPtr4);
		variantsListPtr->getNextCompoundVariant(variantPtr5);
		*/

		ASSERT_EQ(variantPtr1->getPosition(), 249240051);
		ASSERT_STREQ(variantPtr1->getRef().c_str(),    "TG");
		ASSERT_TRUE(std::find(variantPtr1->getAlt().begin(), variantPtr1->getAlt().end(), "TAG") != variantPtr1->getAlt().end());
		ASSERT_EQ(variantPtr1->getAlt().size(), 1);

		/*
		ASSERT_EQ(variantPtr2->getPosition(), 249240099);
		ASSERT_STREQ(variantPtr2->getRef().c_str(),    "TG");
		ASSERT_TRUE(std::find(variantPtr2->getAlt().begin(), variantPtr2->getAlt().end(), "TAG") != variantPtr2->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr2->getAlt().begin(), variantPtr2->getAlt().end(), "TAGG") != variantPtr2->getAlt().end());
		ASSERT_EQ(variantPtr2->getAlt().size(), 2);

		ASSERT_EQ(variantPtr3->getPosition(), 249240303);
		ASSERT_STREQ(variantPtr3->getRef().c_str(),    "TG");
		ASSERT_TRUE(std::find(variantPtr3->getAlt().begin(), variantPtr3->getAlt().end(), "TAG") != variantPtr3->getAlt().end());
		ASSERT_EQ(variantPtr3->getAlt().size(), 1);

		ASSERT_EQ(variantPtr4->getPosition(), 249240537);
		ASSERT_STREQ(variantPtr4->getRef().c_str(),    "GGTG");
		ASSERT_TRUE(std::find(variantPtr4->getAlt().begin(), variantPtr4->getAlt().end(), "GG") != variantPtr4->getAlt().end());
		ASSERT_TRUE(std::find(variantPtr4->getAlt().begin(), variantPtr4->getAlt().end(), "GGGG") != variantPtr4->getAlt().end());
		ASSERT_EQ(variantPtr4->getAlt().size(), 2);
		*/

		// ASSERT_EQ(variantPtr2, nullptr);
	}
}
