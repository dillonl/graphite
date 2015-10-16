#ifndef GRAPHITE_VG_VARIANT_GRAPH_TESTS_HPP
#define GRAPHITE_VG_VARIANT_GRAPH_TESTS_HPP
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
#include "core/variants/VariantListVCFPreloaded.h"
#include "core/reference/FastaReference.h"

 namespace
 {

	 class VGTest : public graphite::vg::VariantGraph
	 {
	 public:
		 typedef std::shared_ptr< VGTest > VGTestPtr;

		 VGTest(graphite::IReference::SharedPtr referencePtr, graphite::IVariantList::SharedPtr variantListPtr) :
			 VariantGraph(referencePtr, variantListPtr)
		 {

		 }

		 ~VGTest()
		 {
		 }

		 bool testGetNextCompoundVariant(graphite::Variant::SharedPtr& variant)
		 {
			 return graphite::IGraph::getNextCompoundVariant(variant);
		 }

	 };

	 //std::string m_reference_string = "GGCCTCAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTGAGATTACAGGTGTGAACCACTGCACCTGGTCTCAAATCAGAAAATCTTTGAATTCATCTAATGGTCACCTATCCTGTGGGCCCTCACTTTGAGATATTTTGCCTTTTTTGGCCAAACCAATATGTAGCCTCCATGTATTATATGACCTTGCCTGCAACCTCTGCCTTCCCACCTTTAAAAACCCTTACACATAAGCCATCAGGGAGATTAGGCCTTAAGGATTAGCTGCCTGATACTCCTTGCTTGCTGCCTGCAATAAATTCCTCAACTTCTGTCTCAGCAATGCCGATATCAGTGCTTGACTTTGATAGGCTGGGTGGGTGGACCCAAATTTGGTTTGGTGACCCTTTGAGCTTAGATTCAAAATTCTAGTTTTGTCACTCTGCAGTTTTGTGATCTTGAGCAAGTTACTTAACCTCTCTGAGCCTTGTTTGTCATGTGTAATGAAAAGAGCTATACTTACCTTGTGAGGTAGTCCTCAGGATTCAATGAGATAATAAGTACTCACTAAACAAAACTCGTTATTACAAAAGAATCACTTTGTCTCTGAAGTGGGCAATTCAACCCATTTCTAGGAGATTTTAAACATGATTTTAGATATTTGGTGTGATTTTGTGAATGGGTTTATCGTTAATAGCTTTCATGCTCCAGAATTTTCTTGAATAATAGGTTTTTGCAAAGTGCATTCCATGGAATACTCATTTGGGTGACGTTAATAGACATCACTCAAAAGCTGGGTGAATATTACAATGTTTACTTCATCTGTAACAAGCTGAGTAGCTACAGTACATATCTAAGAGGGGGCTCTAATTCTCAATATTTTCCAAATTTATTAGATCACAGACTTTTCTTTTAGTGAAGTGCTTAATGAAACTTAAGTTCTGTGAAAAGTACTTTGAGAAATATTGCTTTAAAAAGAAAAAGATTGAGCCCTGTATCAGGGGAAATATCTAATATTATATTAAACAAAAAAGTCCCA";
	 const std::string m_reference_string = "AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG";


	 TEST(VariantGraphTest, getCompoundVariantEmptyTest)
	 {
		 graphite::position pos = 6000;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);
		 ASSERT_TRUE(variantPtr == NULL);
	 }


	 TEST(VariantGraphTest, getCompoundVariantSingleTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 3;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, {"T"});

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		 graphite::Variant::SharedPtr variantPtr;
         variantGraph->testGetNextCompoundVariant(variantPtr);
		 ASSERT_STREQ(variantPtr->getRef().c_str(), std::string(m_reference_string.c_str() + variantPositionOffset, 1).c_str());

         variantGraph->testGetNextCompoundVariant(variantPtr);
		 ASSERT_TRUE(variantPtr == NULL);
	 }

	 TEST(VariantGraphTest, getCompoundVariantSNPSinglePositionTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 3;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);

		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 1, {"A"});
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});


		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);
		 std::string reference = std::string(m_reference_string.c_str() + variantPositionOffset, 1);

		 ASSERT_STREQ(variantPtr->getRef().c_str(), reference.c_str());
		 ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "A");
		 ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "G");
	 }

	 TEST(VariantGraphTest, getCompoundVariantMultiSinglePositionTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 3;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 2, {"AC"});
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);

		 std::string reference = std::string(m_reference_string.c_str() + variantPositionOffset, 2);
		 ASSERT_STREQ(variantPtr->getRef().c_str(), reference.c_str());
		 ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "AC");
		 ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "GA");
	 }

	 TEST(VariantGraphTest, getCompoundVariantMultiWithDeletionTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 3;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 // ref = TAAAAGTTAT
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "1", 10, {"TT"});
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset + 1, "2", 1, {"G"});

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);

		 std::string reference = std::string(m_reference_string.c_str() + variantPositionOffset, 10);
		 ASSERT_STREQ(variantPtr->getRef().c_str(), reference.c_str());
		 ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "TT");
		 // ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "AGGT");
		 ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "TGAAAGTTAT");
	 }

	 TEST(VariantGraphTest, getCompoundVariantMultipleSingleVariantTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("G");
		 variants.push_back("GA");
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, variants);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);
		 std::string reference = std::string(m_reference_string.c_str() + variantPositionOffset, 1);
		 ASSERT_STREQ(variantPtr->getRef().c_str(), reference.c_str());
		 ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "G");
		 ASSERT_STREQ(variantPtr->getAlt()[1].c_str(), "GA");
		 ASSERT_EQ(variantPtr->getAlt().size(), 2);
	 }

	 TEST(VariantGraphTest, getCompoundVariantNoOverlapTest)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos, "1", 1, {"T"});
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, "2", 1, {"G"});

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 graphite::Variant::SharedPtr variantPtr;
		 variantGraph->testGetNextCompoundVariant(variantPtr);
		 std::string reference = std::string(m_reference_string.c_str() + variantPositionOffset, 1);
		 ASSERT_STREQ(variantPtr->getRef().c_str(), reference.c_str());
		 ASSERT_STREQ(variantPtr->getAlt()[0].c_str(), "T");
		 ASSERT_EQ(variantPtr->getAlt().size(), 1);
	 }

	 TEST(VariantGraphTest, ConstructSingleVertexGraph)
	 {
		 graphite::position pos = 6000;
		graphite::position variantPositionOffset = 1;
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
        //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG

		auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		variantGraph->constructGraph();
		// variantGraph->printGraph("test1.dot");
	}

	 TEST(VariantGraphTest, ConstructSingleVariantVertexGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, { "G" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructSinglePositionMultiVariantVertexGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("G");
		 variants.push_back("GA");
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, variants);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructSingleVariantVertexStartPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos, ".", 1, { "G" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructSinglePositionMultiVariantVertexStartPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("G");
		 variants.push_back("GA");
		 testReferenceVariantGenerator.addVariant(pos, ".", 1, variants);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructSinglePositionMultiVariant_MSVVertexStartPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 5;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 10, { "A" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructMultiPositionMultiVariant_MSVVertexStartPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = 5;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("G");
		 variants.push_back("GA");
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 10, variants);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructMultiplePositionMultiVariantVertexGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position firstVariantPositionOffset = 10;
		 graphite::position secondVariantPositionOffset = 19;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("G");
		 variants.push_back("GA");
		 testReferenceVariantGenerator.addVariant(pos + firstVariantPositionOffset, ".", 1, variants);
		 testReferenceVariantGenerator.addVariant(pos + secondVariantPositionOffset, ".", 1, { "G" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructMultiplePositionMultiVariantVertexBackToBackGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position firstVariantPositionOffset = 10;
		 graphite::position secondVariantPositionOffset = 19;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + firstVariantPositionOffset, ".", 1, { "A"});
		 testReferenceVariantGenerator.addVariant(pos + firstVariantPositionOffset + 1, ".", 1, { "G" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructSingleVariantVertexEndPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = m_reference_string.size() - 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, { "A" });

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }

	 TEST(VariantGraphTest, ConstructMultiVariantVertexEndPositionGraph)
	 {
		 graphite::position pos = 6000;
		 graphite::position variantPositionOffset = m_reference_string.size() - 1;
		 graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", pos);
		 //AAATAAAAGTTATCCACCCACCTTGGCCTCCCAAAGCGCTG
		 std::vector< std::string > variants;
		 variants.push_back("T");
		 variants.push_back("A");
		 testReferenceVariantGenerator.addVariant(pos + variantPositionOffset, ".", 1, variants);

		 auto variantGraph = std::make_shared< VGTest >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());
		 variantGraph->constructGraph();
		 // variantGraph->printGraph("test1.dot");
	 }


	TEST(VariantGraphTest, ConstructGraphYChr)
	{

		graphite::Region::SharedPtr regionPtr = std::make_shared< graphite::Region >("20:60343-61343");
		// graphite::Region::SharedPtr regionPtr = std::make_shared< graphite::Region >("20");
		// problem is at position: 2820410
		// graphite::Region::SharedPtr regionPtr = std::make_shared< graphite::Region >("Y:2820000-2821510");

		std::string fastaPath = TEST_FASTA_FILE;
		// std::string vcfPath = TEST_1KG_CHRY_VCF_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;
		// std::string vcfPath = TEST_1KG_VCF_FILE;

		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);
        auto vcfFileReader = std::make_shared< graphite::VariantListVCFPreloaded >(vcfPath, regionPtr);
        auto variantGraph = std::make_shared< VGTest >(fastaReferencePtr, vcfFileReader);
		vcfFileReader->loadVariantsFromFile();
		variantGraph->constructGraph();
		std::cout << "graph constructed" << std::endl;

		variantGraph->printGraph("test1.dot");


		// graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< graphite::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}


}

#endif //VARIANT_GRAPH_TESTS_HPP
