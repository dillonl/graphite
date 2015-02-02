#ifndef GWIZ_VG_VARIANT_GRAPH_TESTS_HPP
#define GWIZ_VG_VARIANT_GRAPH_TESTS_HPP
#include <chrono>
#include <thread>

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

	class VGTest : public gwiz::vg::VariantGraph
	{
	public:
		typedef std::shared_ptr< VGTest > VGTestPtr;

		VGTest(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr)
		{
			m_reference_ptr = referencePtr;
			m_variant_list_ptr = variantListPtr;
			m_graph_ptr = std::make_shared< gwiz::vg::VariantGraph::Graph >();
		}

		~VGTest()
		{
		}

		void constructGraph() override
		{
			constructGraph();
		}

		bool getNextCompoundVariant(gwiz::Variant::SharedPtr& variant) override
		{
			return gwiz::vg::VariantGraph::getNextCompoundVariant(variant);
		}

		gwiz::Variant::SharedPtr buildCompoundVariant(const std::string& referenceString, const std::vector< gwiz::Variant::SharedPtr >& variants) override
		{
			return gwiz::vg::VariantGraph::buildCompoundVariant(referenceString, variants);
		}

		bool doTest(gwiz::Variant::SharedPtr& variant)
		{
			using namespace gwiz;
			std::cout << "dotest " << std::endl;
			// the first time we call this function we need to get the first variant
			if (!m_next_variant_init)
			{
				m_next_variant_init = true;
				m_variant_list_ptr->getNextVariant(this->m_next_variant);
			}
			variant = this->m_next_variant; // set the variant to be returned
			std::cout << "var: " << variant.get() << std::endl;
			if (this->m_next_variant == NULL) { return false; } // if the variant is NULL then return false, this indicates we are past the region of interest

			std::string referenceString = variant->getRef()[0];
			std::vector< Variant::SharedPtr > variants;
			Variant::SharedPtr nextVariant;
			bool variantAdded = false;


			// loop through variants until the variants stop overlapping.
			// As we loop through build a concatenated reference string
			// that represents the entire overlapped variants reference.
			// for those overlapped variants add them to a vector so
			// a "compound variant" can be generated.
			while(m_variant_list_ptr->getNextVariant(nextVariant))
			{
				position variantEndPosition = (variant->getPosition() + referenceString.size());
				if (variantEndPosition < nextVariant->getPosition())
				{
					break;
				}
				// this is a minor efficiency, even though this is a bit ugly
				// it is more efficient to add the variant after checking that this
				// is a compound variant because it is so rare
				if (!variantAdded)
				{
					variants.push_back(variant); // we will build a compound variant with all these variants
					variantAdded = true;
				}
				std::string nextReferenceString = nextVariant->getRef()[0];
				position nextVariantEndPosition = (nextVariant->getPosition() + nextReferenceString.size());
				// if the next variant has reference at a further position then add it to the end of the referenceString
				if (nextVariantEndPosition > variantEndPosition)
				{
					position referenceDelta = (nextVariantEndPosition - variantEndPosition);
					referenceString += nextReferenceString.substr(referenceDelta);
				}
				variants.push_back(nextVariant); // we will build a compound variant with all these variants
			}
			std::cout << "empty: " << variants.size() << std::endl;
			if (!variants.empty())
			{
				variant = buildCompoundVariant(referenceString, variants);
			}
			this->m_next_variant = nextVariant; // set the next variant
			return true;
		}
	};

	//std::string m_reference_string = "GGCCTCAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTGAGATTACAGGTGTGAACCACTGCACCTGGTCTCAAATCAGAAAATCTTTGAATTCATCTAATGGTCACCTATCCTGTGGGCCCTCACTTTGAGATATTTTGCCTTTTTTGGCCAAACCAATATGTAGCCTCCATGTATTATATGACCTTGCCTGCAACCTCTGCCTTCCCACCTTTAAAAACCCTTACACATAAGCCATCAGGGAGATTAGGCCTTAAGGATTAGCTGCCTGATACTCCTTGCTTGCTGCCTGCAATAAATTCCTCAACTTCTGTCTCAGCAATGCCGATATCAGTGCTTGACTTTGATAGGCTGGGTGGGTGGACCCAAATTTGGTTTGGTGACCCTTTGAGCTTAGATTCAAAATTCTAGTTTTGTCACTCTGCAGTTTTGTGATCTTGAGCAAGTTACTTAACCTCTCTGAGCCTTGTTTGTCATGTGTAATGAAAAGAGCTATACTTACCTTGTGAGGTAGTCCTCAGGATTCAATGAGATAATAAGTACTCACTAAACAAAACTCGTTATTACAAAAGAATCACTTTGTCTCTGAAGTGGGCAATTCAACCCATTTCTAGGAGATTTTAAACATGATTTTAGATATTTGGTGTGATTTTGTGAATGGGTTTATCGTTAATAGCTTTCATGCTCCAGAATTTTCTTGAATAATAGGTTTTTGCAAAGTGCATTCCATGGAATACTCATTTGGGTGACGTTAATAGACATCACTCAAAAGCTGGGTGAATATTACAATGTTTACTTCATCTGTAACAAGCTGAGTAGCTACAGTACATATCTAAGAGGGGGCTCTAATTCTCAATATTTTCCAAATTTATTAGATCACAGACTTTTCTTTTAGTGAAGTGCTTAATGAAACTTAAGTTCTGTGAAAAGTACTTTGAGAAATATTGCTTTAAAAAGAAAAAGATTGAGCCCTGTATCAGGGGAAATATCTAATATTATATTAAACAAAAAAGTCCCA";
	const std::string m_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";

	TEST(VariantGraphTest, getCompoundVariantEmptyTest)
	{
		gwiz::position pos = 6000;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		gwiz::Variant::SharedPtr variantPtr;
		variantGraph->getNextCompoundVariant(variantPtr);
		ASSERT_TRUE(variantPtr == NULL);
	}

	TEST(VariantGraphTest, getCompoundVariantSingleTest)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 3;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, {"T"});

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		gwiz::Variant::SharedPtr variantPtr;
		variantGraph->getNextCompoundVariant(variantPtr);
		ASSERT_STREQ(variantPtr->getRef()[0].c_str(), std::string(m_reference_string.c_str() + variantPositionOffset, 1).c_str());

		variantGraph->getNextCompoundVariant(variantPtr);
		// variantGraph->doTest(variantPtr);
		ASSERT_TRUE(variantPtr == NULL);
	}

	TEST(VariantGraphTest, getCompoundVariantSNPSinglePositionTest)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 3;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 1, {"A"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		gwiz::Variant::SharedPtr variantPtr;
		// variantGraph->doTest(variantPtr);
		variantGraph->getNextCompoundVariant(variantPtr);
		ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "A");
		ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "G");
	}

	TEST(VariantGraphTest, getCompoundVariantSNPSinglePositionTest1)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 3;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 1, {"A"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		gwiz::Variant::SharedPtr variantPtr;
		variantGraph->getNextCompoundVariant(variantPtr);
		// variantGraph->doTest(variantPtr);
		ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "A");
		ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "G");
	}

	TEST(VariantGraphTest, getCompoundVariantSNPSinglePositionTest2)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 3;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 1, {"A"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		gwiz::Variant::SharedPtr variantPtr;
		variantGraph->getNextCompoundVariant(variantPtr);
		// variantGraph->doTest(variantPtr);
		ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "A");
		ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "G");
	}

	TEST(VariantGraphTest, getCompoundVariantMultiSinglePositionTest)
	{
		gwiz::position pos = 6000;
		gwiz::position variantPositionOffset = 3;
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 2, {"AC"});
		testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		// std::this_thread::sleep_for (std::chrono::milliseconds(1));
		// std::this_thread::sleep_for (std::chrono::milliseconds(1000));

		// auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		gwiz::Variant::SharedPtr variantPtr;
		variantGraph->getNextCompoundVariant(variantPtr);
		// variantGraph->doTest(variantPtr);

		// ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "AC");
		// ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "GA");
	}

	/*

	TEST_F(VariantGraphTest, ConstructGraphNoVariants)
	{
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
	}

	TEST_F(VariantGraphTest, MultipleVariantsSequentially)
	{
		gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		testReferenceVariantGenerator.addVariant(6001, ".", 1, {"T"});
		testReferenceVariantGenerator.addVariant(6002, ".", 1, {"G"});
		testReferenceVariantGenerator.addVariant(6003, ".", 1, {"A"});
		testReferenceVariantGenerator.addVariant(6004, ".", 1, {"C"});
		testReferenceVariantGenerator.addVariant(6005, ".", 1, {"G"});
		// testReferenceVariantGenerator.addVariant(6015, ".", 1, {"AAAC","GG"});
		// testReferenceVariantGenerator.addVariant(6020, ".", 1, {"A","CCC"});

		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->printGraph("test.dot");
	}
	*/


	/*
	TEST_F(VariantGraphTest, ConstructGraphYChr)
	{

		// gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20:60343-62965354");
		//gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("Y:2655180-2656128");
		// problem is at position: 2820410
		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("Y:2820400-2820410");

		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHRY_VCF_FILE;
		// std::string vcfPath = TEST_1KG_VCF_FILE;

		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReader = std::make_shared<gwiz::VCFFileReader>(vcfPath, regionPtr);
		auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(fastaReferencePtr, vcfFileReader);
		variantGraph->printGraph("test1.dot");


		// gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}
	*/

}

#endif //VARIANT_GRAPH_TESTS_HPP
