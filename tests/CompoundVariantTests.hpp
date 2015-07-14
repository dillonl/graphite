#ifndef GWIZ_COMPOUNDVARIANTTESTS_HPP
#define GWIZ_COMPOUNDVARIANTTESTS_HPP

#include <vector>

#include "TestConfig.h"

#include "core/sequence/SequenceManager.h"
#include "core/allele/Allele.h"
#include "core/allele/EquivalentAllele.h"
#include "core/parser/VCFParser.hpp"
#include "core/variant/Variant.h"
#include "core/variant/VCFFileReader.h"
#include "core/graph/IGraph.h"

namespace
{
	using namespace gwiz;

	void testCompoundVariants(std::vector< std::string > variantLines, std::vector< IAllele::SharedPtr > refAllelePtr, std::vector< std::vector< IAllele::SharedPtr > > altAlleles)
	{
		VariantParser< const char* > vcfParser;
		std::vector< IVariant::SharedPtr > variantPtrs;
		for (auto& variantLine : variantLines)
		{
			variantPtrs.emplace_back(Variant::BuildVariant(variantLine, vcfParser));
		}
		auto variantsListPtr = std::make_shared< VariantList >(variantPtrs);
		variantsListPtr->sort();
		variantsListPtr->normalizeOverlappingVariants();

		uint32_t variantCounter = 0;
		IVariant::SharedPtr variantPtr;
		while (variantsListPtr->getNextVariant(variantPtr))
		{
			auto refMetaDataPtr = variantPtr->getRefAllelePtr()->getAlleleMetaData();
			auto testRefMetaDataPtr = refAllelePtr[variantCounter]->getAlleleMetaData();
			ASSERT_STREQ(variantPtr->getRefAllelePtr()->getSequence(), refAllelePtr[variantCounter]->getSequence());
			ASSERT_EQ(refMetaDataPtr->getPaddingPrefix(), testRefMetaDataPtr->getPaddingPrefix());
			ASSERT_EQ(refMetaDataPtr->getPaddingSuffix(), testRefMetaDataPtr->getPaddingSuffix());
			for (size_t i = 0; i < altAlleles[variantCounter].size(); ++i)
			{
				auto altAllelePtr = altAlleles[variantCounter][i];
				auto altMetaDataPtr = altAllelePtr->getAlleleMetaData();
				ASSERT_STREQ(variantPtr->getAltAllelePtrs()[i]->getSequence(), altAllelePtr->getSequence());
				ASSERT_EQ(altMetaDataPtr->getPaddingPrefix(), altMetaDataPtr->getPaddingPrefix());
				ASSERT_EQ(altMetaDataPtr->getPaddingSuffix(), altMetaDataPtr->getPaddingSuffix());
			}
			ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), altAlleles[variantCounter].size());
			++variantCounter;
		}
	}

	IAllele::SharedPtr makeAllele(const std::string& seq, std::vector< std::tuple< uint16_t, uint16_t > > padding)
	{
		auto sequencePtr = SequenceManager::Instance()->getSequence(seq);
		IAllele::SharedPtr allelePtr;
		if (padding.size() > 1)
		{
			allelePtr = std::make_shared< EquivalentAllele >(sequencePtr);
			for (auto& pad : padding)
			{
				auto metaDataPtr = std::make_shared< AlleleMetaData >(std::get< 0 >(pad), std::get< 1 >(pad));
			}
		}
		else
		{
			allelePtr = std::make_shared< Allele >(sequencePtr);
			auto metaDataPtr = std::make_shared< AlleleMetaData >(std::get< 0 >(padding[0]), std::get< 1 >(padding[0]));
		}
		return allelePtr;
	}

	TEST(CompoundVariantTests, CompoundVariantSimple)
	{
		std::vector< std::string > variantLines = { "20\t60\t.\tTAAGT\tT\t.\t.\t.", "20\t61\t.\tAA\tAGAAG\t.\t.\t.\t.", "20\t63\t.\tGTTAAC\tGTAG,G\t.\t.\t." };
		std::vector< IAllele::SharedPtr > refAllelePtrList = { makeAllele("TAAGTTAAC", {std::make_tuple< uint16_t, uint16_t >(0,0)}) };
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		std::vector< IAllele::SharedPtr > altAllelePtrs;
		altAllelePtrs.emplace_back(makeAllele("TTAAC", {std::make_tuple< uint16_t, uint16_t >(0,4)}));
		altAllelePtrs.emplace_back(makeAllele("TAGAAGGTTAAC", {std::make_tuple< uint16_t, uint16_t >(1,6)}));
		altAllelePtrs.emplace_back(makeAllele("TAAGTAG", {std::make_tuple< uint16_t, uint16_t >(3,0)}));
		altAllelePtrs.emplace_back(makeAllele("TAAG", {std::make_tuple< uint16_t, uint16_t >(3,0)}));
		altAllelePtrsList.emplace_back(altAllelePtrs);
		testCompoundVariants(variantLines, refAllelePtrList, altAllelePtrsList);
	}


	/*
	  TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT
	  TAATATATGTAATATATATTATATATGTAATATAATATATGTAAT
	  *********TAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT
	  **********************************************************GT**************************************************

	  TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT
	  *********TAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT
	  *********TAATATATGT****************************************************************************************************

	*/

	TEST(CompoundVariantTests, BuildCompoundVariant)
	{
		std::vector< std::string > variantLines = { "20\t20301046\t.\tTAATATATGTAATATATATTATATATGTAATATAATATATGTAAT\tT\t.\t.\t.", "20\t20301055\t.\tTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT\tT\t.\t.\t.", "20\t20301104\t.\tGT\tG\t.\t.\t." };
		std::vector< IAllele::SharedPtr > refAllelePtrList = { makeAllele("TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", {std::make_tuple< uint16_t, uint16_t >(0,65),std::make_tuple< uint16_t, uint16_t >(9,0),std::make_tuple< uint16_t, uint16_t >(57,50)}) };
		std::vector< std::vector< IAllele::SharedPtr > > altAllelePtrsList = {};
		std::vector< IAllele::SharedPtr > altAllelePtrs;
		altAllelePtrs.emplace_back(makeAllele("TATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", {std::make_tuple< uint16_t, uint16_t >(0,108)}));
		altAllelePtrs.emplace_back(makeAllele("TAATATATGT", {std::make_tuple< uint16_t, uint16_t >(9,90)}));
		altAllelePtrs.emplace_back(makeAllele("TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT", {std::make_tuple< uint16_t, uint16_t >(3,0)}));
		altAllelePtrsList.emplace_back(altAllelePtrs);
		testCompoundVariants(variantLines, refAllelePtrList, altAllelePtrsList);
		/*
		  std::vector< std::string > refSequence = { "TAATATATGTAATATATATTATATATGTAATATAATATATGTAATATATATTATATATGTAATATATAATATATGTAATATATAATATATGTAATATATATTATATATGT" };
		  std::vector< std::vector< std::string > > altSequences ={ { "", "", "" } };
		std::vector< std::vector< std::tuple< uint32_t, uint32_t > > > altSequencePadding = { { std::make_tuple< uint32_t, uint32_t >(0,65), std::make_tuple< uint32_t, uint32_t >(9,0), std::make_tuple< uint32_t, uint32_t >(58,50) } };
		testCompoundVariants(variantLines, refSequence, altSequences, altSequencePadding);
		*/
	}
/*

	TEST(CompoundVariantTests, BuildSemanticCompoundVariant)
	{
		std::vector< std::string > variantLines = { "1\t105774\t.\tGT\tG\t108.73\ttAC=1;AF=0.500;AN=2;BaseQRankSum=-1.532;ClippingRankSum=0.143;DP=11;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=21.91;MQ0=0;MQRankSum=0.466;QD=9.88;ReadPosRankSum=0.143\tGT:AD:DP:GQ:PL\t0/1:7,5:12:99:167,0,249\t", "1\t105774\t.\tGTT\tGT\t3.22758\t.\tAB=0.277778;ABP=10.7311;AC=1;AF=0.5;AN=2;AO=5;CIGAR=1M1D1M;DP=18;DPB=16.3333;DPRA=0;EPP=6.91895;EPPR=3.17734;GTI=0;LEN=1;MEANALT=1;MQM=21.8;MQMR=16.5385;NS=1;NUMALT=1;ODDS=0.0976763;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=185;QR=465;RO=13;RPP=6.91895;RPPR=3.17734;RUN=1;SAF=5;SAP=13.8677;SAR=0;SRF=13" };
		std::vector< std::string > refSequence = { "GTT" };
		std::vector< std::vector< std::string > > altSequences ={ { "GT" } };
		std::vector< std::vector< std::tuple< uint32_t, uint32_t > > > altSequencePadding = { { std::make_tuple< uint32_t, uint32_t >(0,1) } };
		testCompoundVariants(variantLines, refSequence, altSequences, altSequencePadding);
	}

	TEST(CompoundVariantTests, BuildCompoundVariant3)
	{
		std::vector< std::string > variantLines = { "20\t69506\t.\tG\tGACAC\t590.52\t.\tAC=2;AF=1.00;AN=2;DP=41;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=3.60\tGT:AD:DP:GQ:PL\t1/1:0,25:25:82:1163,82,0", "20\t69506\t.\tGACACACACACACACACACACACACACACA\tGACACACACACACACACACACACACACACACACA\t506.892\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=22;CIGAR=1M4I29M;DP=30;DPB=46.5333;DP" };
		std::vector< std::string > refSequence = { "GACACACACACACACACACACACACACACA" };
		std::vector< std::vector< std::string > > altSequences ={ { "GACACACACACACACACACACACACACACACACA" } };
		std::vector< std::vector< std::tuple< uint32_t, uint32_t > > > altSequencePadding = { { std::make_tuple< uint32_t, uint32_t >(0,29) } };
		testCompoundVariants(variantLines, refSequence, altSequences, altSequencePadding);
	}

	TEST(CompoundVariantTests, BuildGraph4)
	{
		std::vector< std::string > variantLines = { "20\t72104\t.\tTA\tT\t1205.73\t.\tAC=2;AF=1.00;AN=2;DP=39;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.92\tGT:AD:DP:GQ:PL  1/1:0,36:36:99:1243,109,0", "20\t72104\t.\tTAA\tTA\t1062.05\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=36;CIGAR=1M1D1M;DP=39;DPB=59.6667;DPRA=0;EPP=9.04217;EPPR=5.18177;GTI=0;LEN=1;MEANALT=3;MQM=6", "20\t72719\t.\tC\tT\t1779.36\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.63072;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N" };
		std::vector< std::string > refSequence = { "TAA", "C" };
		std::vector< std::vector< std::string > > altSequences ={ { "TA" }, { "T" } };
		std::vector< std::vector< std::tuple< uint32_t, uint32_t > > > altSequencePadding = { { std::make_tuple< uint32_t, uint32_t >(0,1) }, { std::make_tuple< uint32_t, uint32_t >(0,0) } };
		testCompoundVariants(variantLines, refSequence, altSequences, altSequencePadding);
	}

	TEST(CompoundVariantTests, BuildCompoundVariant5)
	{
		std::vector< std::string > variantLines = { "20\t86005\t.\tC\tCA\t31.73\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=-0.361;ClippingRankSum=0.361;DP=41;FS=7.782;MLEAC=1;MLEAF=0.500;MQ=58.77;MQ0=0;MQRankSum=0", "20\t86164\t.\tAG\tA\t2371.73\t.\tAC=2;AF=1.00;AN=2;DP=65;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.22     GT:AD:DP:GQ:PL  1/1:0,65:65:99:2409,196,0", "20\t86164\t.\tAGG\tAG\t2037.92\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=65;CIGAR=1M1D1M;DP=66;DPB=45;DPRA=0;EPP=3.04371;EPPR=0;GTI=0;LEN=1;MEANALT=2;MQM=59.8615;MQMR", "20\t87263\t.\tGA\tG\t1556.73\t.\tAC=2;AF=1.00;AN=2;DP=60;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=25.95     GT:AD:DP:GQ:PL  1/1:0,53:53:99:1594,160,0", "20\t87263\t.\tGAAAAAAAAAT\tGAAAAAAAAT\t1679.52\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=54;CIGAR=1M1D9M;DP=55;DPB=56.4545;DPRA=0;EPP=3.6537;EPPR=0;GTI=0;LEN=1;MEANAL", "20\t87416\t.\tA\tC\t1752.74\t.\tAB=0;ABP=0;AC=2;AF=1;AN=2;AO=56;CIGAR=1X;DP=56;DPB=56;DPRA=0;EPP=3.16541;EPPR=0;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=0;NS=1;N" };
		std::vector< std::string > refSequence = { "C", "AGG", "GAAAAAAAAAT", "A" };
		std::vector< std::vector< std::string > > altSequences ={ { "CA" }, { "AG" }, { "GAAAAAAAAT" }, { "C" } };
		std::vector< std::vector< std::tuple< uint32_t, uint32_t > > > altSequencePadding = { { std::make_tuple< uint32_t, uint32_t >(0,0) }, { std::make_tuple< uint32_t, uint32_t >(0,1) }, { std::make_tuple< uint32_t, uint32_t >(0,9) }, { std::make_tuple< uint32_t, uint32_t >(0,0) } };
		testCompoundVariants(variantLines, refSequence, altSequences, altSequencePadding);
	}
	*/
}

#endif //GWIZ_COMPOUNDVARIANTTESTS_HPP
