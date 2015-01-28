#ifndef GWIZ_VG_VARIANT_GRAPH_TESTS_HPP
#define GWIZ_VG_VARIANT_GRAPH_TESTS_HPP

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "plugins/vg/graph/VariantGraph.h"

#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"

namespace
{


	class VariantGraphTest : public ::testing::Test
	{
	protected:

		//std::string m_reference_string = "GGCCTCAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTGAGATTACAGGTGTGAACCACTGCACCTGGTCTCAAATCAGAAAATCTTTGAATTCATCTAATGGTCACCTATCCTGTGGGCCCTCACTTTGAGATATTTTGCCTTTTTTGGCCAAACCAATATGTAGCCTCCATGTATTATATGACCTTGCCTGCAACCTCTGCCTTCCCACCTTTAAAAACCCTTACACATAAGCCATCAGGGAGATTAGGCCTTAAGGATTAGCTGCCTGATACTCCTTGCTTGCTGCCTGCAATAAATTCCTCAACTTCTGTCTCAGCAATGCCGATATCAGTGCTTGACTTTGATAGGCTGGGTGGGTGGACCCAAATTTGGTTTGGTGACCCTTTGAGCTTAGATTCAAAATTCTAGTTTTGTCACTCTGCAGTTTTGTGATCTTGAGCAAGTTACTTAACCTCTCTGAGCCTTGTTTGTCATGTGTAATGAAAAGAGCTATACTTACCTTGTGAGGTAGTCCTCAGGATTCAATGAGATAATAAGTACTCACTAAACAAAACTCGTTATTACAAAAGAATCACTTTGTCTCTGAAGTGGGCAATTCAACCCATTTCTAGGAGATTTTAAACATGATTTTAGATATTTGGTGTGATTTTGTGAATGGGTTTATCGTTAATAGCTTTCATGCTCCAGAATTTTCTTGAATAATAGGTTTTTGCAAAGTGCATTCCATGGAATACTCATTTGGGTGACGTTAATAGACATCACTCAAAAGCTGGGTGAATATTACAATGTTTACTTCATCTGTAACAAGCTGAGTAGCTACAGTACATATCTAAGAGGGGGCTCTAATTCTCAATATTTTCCAAATTTATTAGATCACAGACTTTTCTTTTAGTGAAGTGCTTAATGAAACTTAAGTTCTGTGAAAAGTACTTTGAGAAATATTGCTTTAAAAAGAAAAAGATTGAGCCCTGTATCAGGGGAAATATCTAATATTATATTAAACAAAAAAGTCCCA";
		std::string m_reference_string = "GGCCTCAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

		// You can remove any or all of the following functions if its body
		// is empty.

		VariantGraphTest()
		{
				// You can do set-up work for each test here.
		}

		virtual ~VariantGraphTest()
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



	TEST_F(VariantGraphTest, ConstructGraph)
	{
		std::string regionString = "20:60808-62964872";
		// std::string regionString = "Y:2655180-2657239";
		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;
		// std::string vcfPath = TEST_1KG_VCF_FILE;
		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >(regionString);

		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(vcfPath);
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(fastaReferencePtr, vcfFileReader);

/*
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);

		testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});

		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
*/
	}
/*
	TEST_F(VariantGraphTest, ConstructGraph)
	{

		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20:60343-62965354");
		//gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("Y:2655180-2657239");
		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;
		// std::string vcfPath = TEST_1KG_VCF_FILE;


		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}
*/
}

#endif //VARIANT_GRAPH_TESTS_HPP
