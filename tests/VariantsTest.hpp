#include <vector>

#include "TestConfig.h"

#include "core/sequence/SequenceManager.h"
#include "core/allele/Allele.h"
#include "core/parser/VCFParser.hpp"
#include "core/variant/Variant.h"
#include "core/variant/VCFFileReader.h"
#include "core/graph/IGraph.h"


namespace
{

	class VariantTest : public graphite::Variant
	{
	public:
		std::string getGenotypeTest() { return getGenotype(); }
		void setAlleleCounts(std::vector< std::string > alleles, std::vector< std::tuple< uint32_t, uint32_t > >& alleleCounts)
		{
			this->m_alt_allele_ptrs.clear();
			this->m_allele_count.clear();
			if (alleles.size() != alleleCounts.size()) { throw "alleles must be matched by allele counts"; }
			for (uint32_t i = 0; i < alleles.size(); ++i)
			{
				auto sequencePtr = graphite::SequenceManager::Instance()->getSequence(alleles[i].c_str());
				auto allelePtr = std::make_shared< graphite::Allele >(sequencePtr);
				this->m_total_allele_count += std::get< 0 >(alleleCounts[i]) + std::get< 1 >(alleleCounts[i]);
				this->m_allele_count[alleles[i]] = alleleCounts[i];
				if (i == 0)	{ m_ref_allele_ptr = allelePtr; }
				else { m_alt_allele_ptrs.emplace_back(allelePtr); }
			}
			this->m_all_allele_ptrs.reserve(this->m_alt_allele_ptrs.size() + 1);
			this->m_all_allele_ptrs.emplace_back(this->m_ref_allele_ptr);
			this->m_all_allele_ptrs.insert(this->m_all_allele_ptrs.end(), this->m_alt_allele_ptrs.begin(), this->m_alt_allele_ptrs.end());
		}
	};

    static std::string VCF_LINE_1 = "Y\t2655180\trs11575897\tG\tA\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	static std::string VCF_LINE_2 = "Y\t2655180\trs11575897\tG\tA,TTA\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line
	static std::string VCF_LINE_3 = "20\t249590\tBI_GS_DEL1_B5_P2733_211\tC\t<CN0>\t100\tPASS\tSVTYPE=DEL;CIEND=-7,7;CIPOS=-7,7;END=250420;CS=DEL_union;MC=EM_DL_DEL10608605;AC=1;AF=0.00019968;NS=2504;AN=5008;EAS_AF=0.0;EUR_AF=0.0;AFR_AF=0.0;AMR_AF=0.0;SAS_AF=0.001\tGT\t0|0"; // is not the complete first line
	static std::string VCF_LINE_4 = "20\t254691\t.\tG\tA\t100\tPASS\tAC=4;AF=0.000798722;AN=5008;NS=2504;DP=14874;EAS_AF=0;AMR_AF=0;AFR_AF=0.003;EUR_AF=0;SAS_AF=0;AA=G|||\tGT\t0|0\t0|0\t0|0\t0|0\t0|0";

	TEST(VariantsTest, ParseVariantChromTest)
	{
		std::string chromVCF = "Y"; // this matches the first variant line of the test_vcf_file
		std::string notChromVCF = "0";

        graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);
		std::string chrom = variantPtr->getChrom();
		EXPECT_STREQ(chrom.c_str(), chromVCF.c_str());
		EXPECT_STRNE(chromVCF.c_str(), notChromVCF.c_str()); // make sure the chrom number and the not chrom number are not equal
		EXPECT_STRNE(chrom.c_str(), notChromVCF.c_str());
	}


	TEST(VariantsTest, ParseVariantPositionTest)
	{
		uint32_t positionVCF = 2655180; // this matches the first variant line of the test_vcf_file
		uint32_t notPositionVCF = 0;

		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

	    uint32_t position = variantPtr->getPosition();
		EXPECT_EQ(position, positionVCF);
		EXPECT_NE(position, notPositionVCF);
		EXPECT_NE(positionVCF, notPositionVCF);
	}

	TEST(VariantsTest, ParseVariantIDTest)
	{
		std::string idVCF = "rs11575897"; // this matches the first variant line of the test_vcf_file
		std::string notIDVCF = "0";

		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

		std::string id = variantPtr->getID();
		EXPECT_STREQ(id.c_str(), idVCF.c_str());
		EXPECT_STRNE(id.c_str(), notIDVCF.c_str()); // make sure the chrom number and the not chrom number are not equal
		EXPECT_STRNE(idVCF.c_str(), notIDVCF.c_str());
	}


	TEST(VariantsTest, ParseVariantRefTest)
	{
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

		auto refAllele = variantPtr->getRefAllelePtr();
		ASSERT_STREQ(refAllele->getSequence(),"G");
		ASSERT_STRNE(refAllele->getSequence(),"A");
	}

	TEST(VariantsTest, ParseVariantAltTest)
	{
		const char* altVCF = "A"; // this matches the first variant line of the test_vcf_file

		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_EQ(altAllelePtrs.size(), 1);
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF);
	}

	TEST(VariantsTest, ParseVariantMultipleAltTest)
	{
		std::vector<std::string> altVCF = {"A","TTA"};

		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_2.c_str(), vcfParser);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF[0].c_str());
		ASSERT_STREQ(altAllelePtrs[1]->getSequence(), altVCF[1].c_str());
	}

	TEST(VariantsTest, ParseVariantMultipleAltDupsTest)
	{
		std::vector<std::string> altVCF = {"A","TTA", "TTA"};

		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_2.c_str(), vcfParser);

		auto altAllelePtrs = variantPtr->getAltAllelePtrs();
		ASSERT_STREQ(altAllelePtrs[0]->getSequence(), altVCF[0].c_str());
		ASSERT_STREQ(altAllelePtrs[1]->getSequence(), altVCF[1].c_str());
		ASSERT_EQ(altAllelePtrs.size(), 2);
	}

	TEST(VariantsTest, ParseVariantQualTest)
	{
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

		std::string qual = variantPtr->getQual();
		ASSERT_STREQ(qual.c_str(),"34439.5");;
	}

	TEST(VariantsTest, ParseVariantFilterTest)
	{
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

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
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_1.c_str(), vcfParser);

		auto infoFields = variantPtr->getInfoFields();
		for (auto infoFieldIter : infoMap)
		{
			ASSERT_TRUE(infoFields.find(infoFieldIter.first) != infoFields.end());
			ASSERT_STREQ(infoFields[infoFieldIter.first].c_str(), infoFieldIter.second.c_str());
		}
		ASSERT_EQ(infoFields.size(), 11);
	}

	TEST(VariantsTest, ParseVariantSymbolidTest)
	{
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_3.c_str(), vcfParser);

		ASSERT_STREQ(variantPtr->getRef().c_str(), "C");
		ASSERT_STREQ(variantPtr->getAltAllelePtrs()[0]->getSequence(), "<CN0>");
	}

	TEST(VariantsTest, ParseVariantQual2Test)
	{
		graphite::VariantParser< const char* > vcfParser;
		graphite::Variant::SharedPtr variantPtr;
		variantPtr = graphite::Variant::BuildVariant(VCF_LINE_4.c_str(), vcfParser);

		std::string qual = variantPtr->getQual();
		ASSERT_STREQ(qual.c_str(),"100");
	}


	/*
	TEST(VariantsTest, TestGetGenotypeSimpleNone)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("./.", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoRefInsufficientCountAlt)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(2,2));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(1,1));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("0/0", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoRefInsufficientPercentAlt)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(50,50));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(29,0));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("0/0", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoRef)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(10,10));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(0,0));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("0/0", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoAlt)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(0,0));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(10,10));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("1/1", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoAltInsufficientCountRef)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(1,1));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(2,2));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("1/1", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoAltInsufficientPercentRef)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(29,0));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(50,50));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("1/1", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeSimpleHet)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(10,10));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(10,10));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("0/1", variantTest.getGenotypeTest().c_str());
	}

	TEST(VariantsTest, TestGetGenotypeHomoSecondAlt)
	{
		VariantTest variantTest;
		std::vector< std::string > alleles;
		std::vector< std::tuple< uint32_t, uint32_t > > alleleCounts;
		alleles.emplace_back("A");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(0,0));
		alleles.emplace_back("AT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(0,0));
		alleles.emplace_back("GAT");
		alleleCounts.emplace_back(std::make_tuple< uint32_t, uint32_t >(50,50));
		variantTest.setAlleleCounts(alleles, alleleCounts);
		ASSERT_STREQ("2/2", variantTest.getGenotypeTest().c_str());
	}
	*/
}
