#ifndef GWIZ_COMPOUNDVARIANTTESTS_HPP
#define GWIZ_COMPOUNDVARIANTTESTS_HPP

#include <vector>

#include "TestConfig.h"

#include "core/sequence/SequenceManager.h"
#include "core/allele/IAllele.h"
#include "core/parser/VCFParser.hpp"
#include "core/variant/Variant.h"
#include "core/variant/VCFFileReader.h"
#include "core/graph/IGraph.h"

namespace
{

	TEST(CompoundVariantTests, CompoundVariantSimple)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t60\t.\tTAAGT\tT\t.\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t61\t.\tAA\tAGAAG\t.\t.\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t63\t.\tGTTAAC\tGTAG,G\t.\t.\t.", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantPtr1);
		variantsVec.emplace_back(variantPtr2);
		variantsVec.emplace_back(variantPtr3);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr;
		variantsListPtr->getNextVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), "TAAGTTAAC");
		std::vector< std::string > altSequences = { "TTAAC", "TAGAAGGTTAAC", "TAAGTAG", "TAAG" };
		for (size_t i = 0; i < altSequences.size(); ++i)
		{
			ASSERT_STREQ(variantPtr->getAltAllelePtrs()[i]->getSequence(), altSequences[i].c_str());
		}
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altSequences.size());
	}

	TEST(CompoundVariantTests, BuildCompoundVariant)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t20301046\t.\tTAATATATGTAATATATATTATATATGTAATATAATATATGTAAT\tT\t.\t.\t.", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t20301055\t.\tTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT\tT\t.\t.\t.", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t20301104\t.\tGT\tG\t.\t.\t.", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantPtr1);
		variantsVec.emplace_back(variantPtr2);
		variantsVec.emplace_back(variantPtr3);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr;
		variantsListPtr->getNextVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT");
		std::vector< std::string > altSequences = { "TATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", "TAATATATGT", "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT" };
		for (size_t i = 0; i < altSequences.size(); ++i)
		{
			ASSERT_STREQ(variantPtr->getAltAllelePtrs()[i]->getSequence(), altSequences[i].c_str());
		}
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altSequences.size());
	}

	TEST(CompoundVariantTests, BuildCompoundVariant2)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("1\t105774\t.\tGT\tG\t108.73\ttAC=1;AF=0.500;AN=2;BaseQRankSum=-1.532;ClippingRankSum=0.143;DP=11;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.91;MQ0=0;MQRankSum=0.466;QD=9.88;ReadPosRankSum=0.143\tGT:AD:DP:GQ:PL\t0/1:7,5:12:99:167,0,249\t", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("1\t105774\t.\tGTT\tGT\t3.22758\t.\tAB=0.277778;ABP=10.7311;AC=1;AF=0.5;AN=2;AO=5;CIGAR=1M1D1M;DP=18;DPB=16.3333;DPRA=0;EPP=6.91895;EPPR=3.17734;GTI=0;LEN=1;MEANALT=1;MQM=21.8;MQMR=16.5385;NS=1;NUMALT=1;ODDS=0.0976763;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=185;QR=465;RO=13;RPP=6.91895;RPPR=3.17734;RUN=1;SAF=5;SAP=13.8677;SAR=0;SRF=13", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantPtr1);
		variantsVec.emplace_back(variantPtr2);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr;
		variantsListPtr->getNextVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), "GTT");
		std::vector< std::string > altSequences = { "GT" };
		for (size_t i = 0; i < altSequences.size(); ++i)
		{
			ASSERT_STREQ(variantPtr->getAltAllelePtrs()[i]->getSequence(), altSequences[i].c_str());
		}
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altSequences.size());
	}

	TEST(CompoundVariantTests, BuildCompoundVariant3)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t69506\t.\tG\tGACAC\t590.52\t.\tAC=2;AF=1.00;AN=2;DP=41;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=3.60\tGT:AD:DP:GQ:PL\t1/1:0,25:25:82:1163,82,0", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t69506\t.\tGACACACACACACACACACACACACACACA\tGACACACACACACACACACACACACACACACACA\t506.892\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=22;CIGAR=1M4I29M;DP=30;DPB=46.5333;DP", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantPtr1);
		variantsVec.emplace_back(variantPtr2);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr;
		variantsListPtr->getNextVariant(variantPtr);

		ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), "GACACACACACACACACACACACACACACA");
		std::vector< std::string > altSequences = { "GACACACACACACACACACACACACACACACACA" };
		for (size_t i = 0; i < altSequences.size(); ++i)
		{
			ASSERT_STREQ(variantPtr->getAltAllelePtrs()[i]->getSequence(), altSequences[i].c_str());
		}
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altSequences.size());
	}

	TEST(CompoundVariantTests, BuildGraph4)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantPtr1 = gwiz::Variant::BuildVariant("20\t72104\t.\tTA\tT\t1205.73\t.\tAC=2;AF=1.00;AN=2;DP=39;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.92\tGT:AD:DP:GQ:PL  1/1:0,36:36:99:1243,109,0", vcfParser);
		auto variantPtr2 = gwiz::Variant::BuildVariant("20\t72104\t.\tTAA\tTA\t1062.05\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=36;CIGAR=1M1D1M;DP=39;DPB=59.6667;DPRA=0;EPP=9.04217;EPPR=5.18177;GTI=0;LEN=1;MEANALT=3;MQM=6", vcfParser);
		auto variantPtr3 = gwiz::Variant::BuildVariant("20\t72719\t.\tC\tT\t1779.36\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.63072;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantPtr1);
		variantsVec.emplace_back(variantPtr2);
		variantsVec.emplace_back(variantPtr3);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr;
		gwiz::IVariant::SharedPtr variantPtrTest;
		gwiz::IVariant::SharedPtr variantPtrNULL;
		variantsListPtr->getNextVariant(variantPtr);
		variantsListPtr->getNextVariant(variantPtrTest);
		variantsListPtr->getNextVariant(variantPtrNULL);

		ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), "TAA");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "TA");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
		ASSERT_EQ(variantPtr->getPosition(), 72104);

		ASSERT_STREQ(variantPtrTest->getRefAllelePtr()->getSequence(),    "C");
		ASSERT_STREQ(variantPtrTest->getAltAllelePtrs()[0]->getSequence(), "T");
		ASSERT_EQ(variantPtrTest->getAltAllelePtrs().size(), 1);
		ASSERT_EQ(variantPtrTest->getPosition(), 72719);

		ASSERT_EQ(variantPtrNULL, nullptr);
	}

	TEST(CompoundVariantTests, BuildCompoundVariant5)
	{
		gwiz::VariantParser< const char* > vcfParser;
		auto variantLinePtr1 = gwiz::Variant::BuildVariant("20\t86005\t.\tC\tCA\t31.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=-0.361;ClippingRankSum=0.361;DP=41;FS=7.782;MLEAC=1;MLEAF=0.500;MQ=58.77;MQ0=0;MQRankSum=0", vcfParser);
		auto variantLinePtr2 = gwiz::Variant::BuildVariant("20\t86164\t.\tAG\tA\t2371.73\t.\tAC=2;AF=1.00;AN=2;DP=65;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.22     GT:AD:DP:GQ:PL  1/1:0,65:65:99:2409,196,0", vcfParser);
		auto variantLinePtr3 = gwiz::Variant::BuildVariant("20\t86164\t.\tAGG\tAG\t2037.92\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=65;CIGAR=1M1D1M;DP=66;DPB=45;DPRA=0;EPP=3.04371;EPPR=0;GTI=0;LEN=1;MEANALT=2;MQM=59.8615;MQMR", vcfParser);
		auto variantLinePtr4 = gwiz::Variant::BuildVariant("20\t87263\t.\tGA\tG\t1556.73\t.\tAC=2;AF=1.00;AN=2;DP=60;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=25.95     GT:AD:DP:GQ:PL  1/1:0,53:53:99:1594,160,0", vcfParser);
		auto variantLinePtr5 = gwiz::Variant::BuildVariant("20\t87263\t.\tGAAAAAAAAAT\tGAAAAAAAAT\t1679.52\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=54;CIGAR=1M1D9M;DP=55;DPB=56.4545;DPRA=0;EPP=3.6537;EPPR=0;GTI=0;LEN=1;MEANAL", vcfParser);
		auto variantLinePtr6 = gwiz::Variant::BuildVariant("20\t87416\t.\tA\tC\t1752.74\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.16541;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N", vcfParser);

		std::vector< gwiz::IVariant::SharedPtr > variantsVec;
		variantsVec.emplace_back(variantLinePtr1);
		variantsVec.emplace_back(variantLinePtr2);
		variantsVec.emplace_back(variantLinePtr3);
		variantsVec.emplace_back(variantLinePtr4);
		variantsVec.emplace_back(variantLinePtr5);
		variantsVec.emplace_back(variantLinePtr6);
		auto variantsListPtr = std::make_shared< gwiz::VariantList >(variantsVec);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		gwiz::IVariant::SharedPtr variantPtr1;
		gwiz::IVariant::SharedPtr variantPtr2;
		gwiz::IVariant::SharedPtr variantPtr3;
		gwiz::IVariant::SharedPtr variantPtr4;
		gwiz::IVariant::SharedPtr variantPtr5;
		variantsListPtr->getNextVariant(variantPtr1);
		variantsListPtr->getNextVariant(variantPtr2);
		variantsListPtr->getNextVariant(variantPtr3);
		variantsListPtr->getNextVariant(variantPtr4);
		variantsListPtr->getNextVariant(variantPtr5);

		ASSERT_EQ(variantPtr1->getPosition(), 86005);
		ASSERT_STREQ(variantPtr1->getRefAllelePtr()->getSequence(),    "C");
		ASSERT_STREQ(variantPtr1->getAltAllelePtrs()[0]->getSequence(), "CA");
		ASSERT_EQ(variantPtr1->getAltAllelePtrs().size(), 1);

		ASSERT_EQ(variantPtr2->getPosition(), 86164);
		ASSERT_STREQ(variantPtr2->getRefAllelePtr()->getSequence(),    "AGG");
		ASSERT_STREQ(variantPtr2->getAltAllelePtrs()[0]->getSequence(), "AG");
		ASSERT_EQ(variantPtr2->getAltAllelePtrs().size(), 1);

		ASSERT_EQ(variantPtr3->getPosition(), 87263);
		ASSERT_STREQ(variantPtr3->getRefAllelePtr()->getSequence(),    "GAAAAAAAAAT");
		ASSERT_STREQ(variantPtr3->getAltAllelePtrs()[0]->getSequence(), "GAAAAAAAAT");
		ASSERT_EQ(variantPtr3->getAltAllelePtrs().size(), 1);

		ASSERT_EQ(variantPtr4->getPosition(), 87416);
		ASSERT_STREQ(variantPtr4->getRefAllelePtr()->getSequence(),    "A");
		ASSERT_STREQ(variantPtr4->getAltAllelePtrs()[0]->getSequence(), "C");
		ASSERT_EQ(variantPtr4->getAltAllelePtrs().size(), 1);

		ASSERT_EQ(variantPtr5, nullptr);
	}
}

#endif //GWIZ_COMPOUNDVARIANTTESTS_HPP
