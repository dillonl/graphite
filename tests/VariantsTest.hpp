#include <vector>

#include "TestConfig.h"

#include "TestClasses.hpp"
#include "core/reference/FastaReference.h"
#include "core/allele/Allele.h"
#include "core/variant/Variant.h"
#include "core/variant/VCFHeader.h"
#include "core/variant/VCFFileReader.h"
#include "core/graph/IGraph.h"


namespace
{
    static std::string VCF_LINE_1 = "Y\t2655180\trs11575897\tG\tA\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	static std::string VCF_LINE_2 = "Y\t2655180\trs11575897\tG\tA,TTA\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	static std::string VCF_LINE_3 = "20\t249590\tBI_GS_DEL1_B5_P2733_211\tC\t<CN0>\t100\tPASS\tSVTYPE=DEL;CIEND=-7,7;CIPOS=-7,7;END=250420;CS=DEL_union;MC=EM_DL_DEL10608605;AC=1;AF=0.00019968;NS=2504;AN=5008;EAS_AF=0.0;EUR_AF=0.0;AFR_AF=0.0;AMR_AF=0.0;SAS_AF=0.001\tGT\t0|0"; // is not the complete first line
	static std::string VCF_LINE_4 = "20\t254691\t.\tG\tA\t100\tPASS\tAC=4;AF=0.000798722;AN=5008;NS=2504;DP=14874;EAS_AF=0;AMR_AF=0;AFR_AF=0.003;EUR_AF=0;SAS_AF=0;AA=G|||\tGT\t0|0\t0|0\t0|0\t0|0\t0|0";
	static std::string VCF_LINE_5 = "1\t75148200\tWG:DEL:a3bcba65\tA\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=-20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47"; // is not the complete first line

	std::string getVariantLineRefOfSize(size_t refSize)
	{
		std::string vcfLine = "1\t75148200\tWG:DEL:a3bcba65\t";

		vcfLine += std::string(refSize, 'T');

		vcfLine += "\tA\t.\t.\tSVTYPE=DEL;SVLEN=-20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47"; // is not the complete first line
		return vcfLine;
	}

	std::string getVariantLineAltOfSize(size_t altSize)
	{
		std::string vcfLine = "1\t75148200\tWG:DEL:a3bcba65\tA\t";

		vcfLine += std::string(altSize, 'T');

		vcfLine += "\t.\t.\tSVTYPE=DEL;SVLEN=-20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47"; // is not the complete first line
		return vcfLine;
	}
/*
	TEST(VariantsTest, TestShouldNotSkipSmallReference)
	{
		uint32_t alleleSizeThreshold = 3000;
		std::string vcfLine = getVariantLineRefOfSize(alleleSizeThreshold);
		graphite::Variant::SharedPtr variantPtr;
        variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), nullptr, 300);
		EXPECT_FALSE(variantPtr->shouldSkip());
	}

	TEST(VariantsTest, TestShouldSkipLargeReference)
	{
		uint32_t alleleSizeThreshold = 3000;
		std::string vcfLine = getVariantLineRefOfSize(alleleSizeThreshold + 1);
		graphite::Variant::SharedPtr variantPtr;
        variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), nullptr, 300);
		EXPECT_TRUE(variantPtr->shouldSkip());
	}

	TEST(VariantsTest, TestShouldSkipLargeAlternate)
	{
		uint32_t alleleSizeThreshold = 3000;
		std::string vcfLine = getVariantLineAltOfSize(alleleSizeThreshold + 1);
		graphite::Variant::SharedPtr variantPtr;
        variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), nullptr, alleleSizeThreshold, 300);
		EXPECT_TRUE(variantPtr->shouldSkip());
	}
*/

	TEST(VariantsTest, ParseVariantChromTest)
	{
		std::string chromVCF = "Y"; // this matches the first variant line of the test_vcf_file
		std::string notChromVCF = "0";

		graphite::Variant::SharedPtr variantPtr;
        variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);
		std::string chrom = variantPtr->getChrom();
		EXPECT_STREQ(chrom.c_str(), chromVCF.c_str());
		EXPECT_STRNE(chromVCF.c_str(), notChromVCF.c_str()); // make sure the chrom number and the not chrom number are not equal
		EXPECT_STRNE(chrom.c_str(), notChromVCF.c_str());
	}


	TEST(VariantsTest, ParseVariantPositionTest)
	{
		uint32_t positionVCF = 2655180; // this matches the first variant line of the test_vcf_file
		uint32_t notPositionVCF = 0;

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

	    uint32_t position = variantPtr->getPosition();
		EXPECT_EQ(position, positionVCF);
		EXPECT_NE(position, notPositionVCF);
		EXPECT_NE(positionVCF, notPositionVCF);
	}

	TEST(VariantsTest, ParseVariantIDTest)
	{
		std::string idVCF = "rs11575897"; // this matches the first variant line of the test_vcf_file
		std::string notIDVCF = "0";

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		std::string id = variantPtr->getID();
		EXPECT_STREQ(id.c_str(), idVCF.c_str());
		EXPECT_STRNE(id.c_str(), notIDVCF.c_str()); // make sure the chrom number and the not chrom number are not equal
		EXPECT_STRNE(idVCF.c_str(), notIDVCF.c_str());
	}


	TEST(VariantsTest, ParseVariantRefTest)
	{
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		auto refAllele = variantPtr->getRefAllelePtr();
		ASSERT_STREQ(refAllele->getSequence(),"G");
		ASSERT_STRNE(refAllele->getSequence(),"A");
	}

	TEST(VariantsTest, ParseVariantAltTest)
	{
		const char* altVCF = "A"; // this matches the first variant line of the test_vcf_file

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_EQ(altAllelePtrs.size(), 1);
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF);
	}

	TEST(VariantsTest, ParseVariantMultipleAltTest)
	{
		std::vector<std::string> altVCF = {"A","TTA"};

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_2.c_str(), nullptr, 300);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF[0].c_str());
		ASSERT_STREQ(altAllelePtrs[1]->getSequence(), altVCF[1].c_str());
	}

	TEST(VariantsTest, ParseSymbolicVariantAltTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", 75148190,75148250, graphite::Region::BASED::ONE);
		std::string sequence = "TGAAGGCCAAAATTCAGATTCAGGACCCCTCCCGGGTAAAAATATATATA";
		auto referencePtr = std::make_shared< graphite::Reference >(sequence, regionPtr);

		// const char* refVCF = "AAATTCAGATTCAGGACCCCT"; // this matches the first variant line of the test_vcf_file
		const char* refVCF = "AATTCAGATTCAGGACCCCTC"; // this matches the first variant line of the test_vcf_file

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_5.c_str(), referencePtr, 300);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		auto refAllele = variantPtr->getRefAllelePtr();
		ASSERT_EQ(altAllelePtrs.size(), 1);
		ASSERT_STREQ(refAllele->getSequence(), refVCF);
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(),"A");
	}

	TEST(VariantsTest, ParseVariantMultipleAltDupsTest)
	{
		std::vector<std::string> altVCF = {"A","TTA", "TTA"};

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_2.c_str(), nullptr, 300);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF[0].c_str());
		ASSERT_STREQ(altAllelePtrs[1]->getSequence(), altVCF[1].c_str());
		ASSERT_EQ(altAllelePtrs.size(), 2);
	}

	TEST(VariantsTest, ParseVariantQualTest)
	{
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		std::string qual = variantPtr->getQual();
		ASSERT_STREQ(qual.c_str(),"34439.5");;
	}

	TEST(VariantsTest, ParseVariantFilterTest)
	{
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		std::string filter = variantPtr->getFilter();
		ASSERT_STREQ(filter.c_str(),"PASS");
	}

	TEST(VariantsTest, ParseVariantInfoTest)
	{
		std::map< std::string, std::string > infoMap;
		infoMap["AA"] = "G";
		infoMap["AC"] = "22";
		infoMap["AF"] = "0.0178427";
		infoMap["AN"] = "1233";
		infoMap["DP"] = "84761";
		infoMap["NS"] = "1233";
		infoMap["AMR_AF"] = "0.0000";
		infoMap["AFR_AF"] = "0.0000";
		infoMap["EUR_AF"] = "0.0000";
		infoMap["SAS_AF"] = "0.0000";
		infoMap["EAS_AF"] = "0.0451";
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), nullptr, 300);

		auto infoFields = variantPtr->getInfoFields();
		for (auto infoFieldIter : infoMap)
		{
			ASSERT_TRUE(infoFields.find(infoFieldIter.first) != infoFields.end());
			ASSERT_STREQ(infoFields[infoFieldIter.first].c_str(), infoFieldIter.second.c_str());
		}
		ASSERT_EQ(infoFields.size(), 11);
	}

	/*
	TEST(VariantsTest, ParseVariantSymbolicTest)
	{
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_3.c_str(), nullptr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "C");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "<CN0>");
	}
	*/

	TEST(VariantsTest, ParseVariantQual2Test)
	{
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_4.c_str(), nullptr, 300);

		std::string qual = variantPtr->getQual();
		ASSERT_STREQ(qual.c_str(),"100");
	}

	static std::string VCF_LINE_SIMPLE = "1\t20\tWG:a3bcba65\tA\tT\t.\t.\tID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;LID=NA12878;RID=NA12878;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47"; // is not the complete first line
	TEST(VariantsTest, ParseVariantSimpleTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_SIMPLE.c_str(), referencePtr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "A");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "T");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	static std::string VCF_LINE_DUP = "1\t20\tWG:DUP:a3bcba65\tA\t<DUP>\t.\t.\tSVTYPE=DUP;SVLEN=20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47"; // is not the complete first line
	TEST(VariantsTest, ParseVariantDupTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_DUP.c_str(), referencePtr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "ACCAAACGTCGTTAGGCCAGT");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "ACCAAACGTCGTTAGGCCAGTCCAAACGTCGTTAGGCCAGT");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	TEST(VariantsTest, ParseVariantDupSVTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_DUP.c_str(), referencePtr, 6);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "ACCAAACNNNNNNGCCAGT");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "ACCAAACNNNNNNGCCAGTCCAAACNNNNNNGCCAGT");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	static std::string VCF_LINE_DEL = "1\t10\tWG:DEL:a3bcba65\tT\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47";
	TEST(VariantsTest, ParseVariantDelTest)
	{
		// CTATGATGTTGATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACA
		// TGATGGAACTGACCAAACGTC
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_DEL.c_str(), referencePtr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "TGATGGAACTGACCAAACGTC");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "T");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	TEST(VariantsTest, ParseVariantDelSVTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_DEL.c_str(), referencePtr, 6);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "TGATGGANNNNNNAACGTC");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "T");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	static std::string VCF_LINE_INV = "1\t10\tWG:INV:a3bcba65\tT\t<INV>\t.\t.\tSVTYPE=INV;SVLEN=20;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47";
	TEST(VariantsTest, ParseVariantInvTest)
	{
		// TGATGGAACTGACCAAACGTC
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_INV.c_str(), referencePtr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "TGATGGAACTGACCAAACGTC");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "TCTGCAAACCAGTCAAGGTAG");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	TEST(VariantsTest, ParseVariantInvSVTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_INV.c_str(), referencePtr, 6);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "TGATGGANNNNNNAACGTC");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "TCTGCAANNNNNNAGGTAG");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	static std::string VCF_LINE_INS = "1\t10\tWG:INS:a3bcba65\tT\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=20;SEQ=AAAAAAAAAAAAAAAAAAAA;ID=a3bcba65;SUPPORT=21,22;MERGED=0;REFINED=1;END=75148220;POS=75148200,75148220;LID=NA12878;RID=NA12878;CIPOS=-10,10;CIEND=-10,10;COLLAPSED=0\tGT:GL:AS:RS\t0/1:-2655.97,-1567.88,-6335.32:24:47";
	TEST(VariantsTest, ParseVariantInsTest)
	{
		// TGATGGAACTGACCAAACGTC
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_INS.c_str(), referencePtr, 300);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "T");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "TAAAAAAAAAAAAAAAAAAAA");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

	TEST(VariantsTest, ParseVariantInsSVTest)
	{
		auto regionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
		auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, regionPtr);

		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_INS.c_str(), referencePtr, 6);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "T");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "TAAAAAANNNNNNAAAAAA");
		ASSERT_EQ(variantPtr->getAltAllelePtrs().size(), 1);
	}

}
