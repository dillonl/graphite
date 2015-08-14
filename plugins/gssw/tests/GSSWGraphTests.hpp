#ifndef GRAPHITE_GSSW_GSSWGRAPH_TESTS_HPP
#define GRAPHITE_GSSW_GSSWGRAPH_TESTS_HPP

#include <chrono>
#include <thread>
#include <memory>

#include "gssw/gssw.h"
#include "gssw/graph/GraphManager.h"

#include "gtest/gtest.h"

#include "config/TestConfig.h"
#include "tests/classes/TestReference.hpp"
#include "tests/classes/TestVariantList.hpp"
#include "tests/classes/TestReferenceVariantGenerator.hpp"
#include "tests/classes/TestAlignmentReader.hpp"
#include "tests/classes/TestAlignmentReaderManager.hpp"
#include "plugins/gssw/graph/GSSWGraph.h"
#include "plugins/gssw/graph/AlignmentReporter.h"
#include "plugins/gssw/graph/GSSWAdjudicator.h"

#include "core/alignment/BamAlignmentReader.h"
#include "core/alignment/BamAlignmentReaderManager.h"
#include "core/alignment/BamAlignmentReaderPreloadManager.h"
#include "core/variant/VCFFileReader.h"
#include "core/variant/IVariant.h"
#include "core/variant/VariantListVCFPreloaded.h"
#include "core/reference/FastaReference.h"
#include "core/util/ThreadPool.hpp"


namespace
{

	class GSSWGraphTest : public graphite::gssw::GSSWGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraphTest > GSSWGraphTestPtr;

		/*
		GSSWGraphTest(graphite::IReference::SharedPtr referencePtr, graphite::IVariantList::SharedPtr variantListPtr) :
			graphite::gssw::GSSWGraph(referencePtr, variantListPtr)
		{
		}

		~GSSWGraphTest()
		{
		}
		*/
		/*
		graphite::gssw::GSSWGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(graphite::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		graphite::gssw::GSSWGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}
		*/

		/*
		static void GenerateGenericGraphTestData(graphite::IReference::SharedPtr& referencePtr, graphite::IVariantList::SharedPtr& variantListPtr)
		{
		}
		*/

		/*
		static void Build31VariantsTestData(graphite::IReference::SharedPtr& referencePtr, graphite::IVariantList::SharedPtr& variantListPtr, graphite::testing::TestAlignmentReader::SharedPtr& alignmentReaderPtr)
		{
			std::string referenceSequence = "ATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAG";

			graphite::position startPosition = 200000;
			std::string referenceString = std::string(startPosition, 'x');
			referenceString += referenceSequence;
			graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceString, "20", 1);
            testReferenceVariantGenerator.addVariant(startPosition + 14, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 18, ".", 1, {"GAAAG"});
			testReferenceVariantGenerator.addVariant(startPosition + 26, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 30, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 104, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 119, ".", 6, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 142, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 220, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 235, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 310, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 315, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 383, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 425, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 431, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 478, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 482, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 492, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 541, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 590, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 652, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 654, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 655, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 658, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 671, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 721, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 755, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 808, ".", 1, {"G"});
			testReferenceVariantGenerator.addVariant(startPosition + 869, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 876, ".", 1, {"A"});
			testReferenceVariantGenerator.addVariant(startPosition + 904, ".", 1, {"T"});
			testReferenceVariantGenerator.addVariant(startPosition + 969, ".", 1, {"C"});
			testReferenceVariantGenerator.addVariant(startPosition + 981, ".", 1, {"C"});

			alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
			alignmentReaderPtr->addAlignment(199924, "GGACCTTTCTCAGCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAGAAAGCTCTAGG");
			alignmentReaderPtr->addAlignment(199933, "TCAGCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGGCAATGCAGC");
			alignmentReaderPtr->addAlignment(199936, "GCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCAACTCTCCCCCAACCCCAAGCAAAGCAGCCTG");
			alignmentReaderPtr->addAlignment(199938, "AGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGG");
			alignmentReaderPtr->addAlignment(199951, "GGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTC");
			alignmentReaderPtr->addAlignment(199970, "AAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGG");
			alignmentReaderPtr->addAlignment(199975, "CAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGGAGAAGAGGGAAGA");
			alignmentReaderPtr->addAlignment(200004, "CACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGG");
			alignmentReaderPtr->addAlignment(200017, "ACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCC");
			alignmentReaderPtr->addAlignment(200042, "AGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTT");
			alignmentReaderPtr->addAlignment(200046, "GACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTA");
			alignmentReaderPtr->addAlignment(200047, "CCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCA");
			alignmentReaderPtr->addAlignment(200047, "CCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCA");
			alignmentReaderPtr->addAlignment(200050, "TCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCA");
			alignmentReaderPtr->addAlignment(200054, "CTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGT");
			alignmentReaderPtr->addAlignment(200066, "AGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGT");
			alignmentReaderPtr->addAlignment(200070, "GGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCC");
			alignmentReaderPtr->addAlignment(200109, "AGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCA");
			alignmentReaderPtr->addAlignment(200163, "TGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAAC");
			alignmentReaderPtr->addAlignment(200183, "CTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGA");
			alignmentReaderPtr->addAlignment(200184, "TTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAA");
			alignmentReaderPtr->addAlignment(200214, "GAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAA");
			alignmentReaderPtr->addAlignment(200215, "AGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAG");
			alignmentReaderPtr->addAlignment(200232, "TTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCGAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACT");
			alignmentReaderPtr->addAlignment(200232, "TTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACT");
			alignmentReaderPtr->addAlignment(200240, "AAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGA");
			alignmentReaderPtr->addAlignment(200252, "AGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTT");
			alignmentReaderPtr->addAlignment(200263, "TGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTC");
			alignmentReaderPtr->addAlignment(200267, "GCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCAC");
			alignmentReaderPtr->addAlignment(200268, "CCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACC");
			alignmentReaderPtr->addAlignment(200292, "AAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAAT");
			alignmentReaderPtr->addAlignment(200345, "CAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCT");
			alignmentReaderPtr->addAlignment(200351, "TTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAA");
			alignmentReaderPtr->addAlignment(200358, "CTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGAT");
			alignmentReaderPtr->addAlignment(200410, "ACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCA");
			alignmentReaderPtr->addAlignment(200410, "ACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCA");
			alignmentReaderPtr->addAlignment(200425, "CTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAA");
			alignmentReaderPtr->addAlignment(200425, "CTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAA");
			alignmentReaderPtr->addAlignment(200434, "CCCTGAAAAGCCTTCCCAAAAAGGATGGGCGCAAACACGATAAAATTGCGACGGCAACGATGTGCAGCTGACTCATCATAACCCGGACATAAACTGCATAC");
			alignmentReaderPtr->addAlignment(200458, "TGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTGAGAACATCCAGGGA");
			alignmentReaderPtr->addAlignment(200484, "AGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAA");
			alignmentReaderPtr->addAlignment(200517, "AGACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGA");
			alignmentReaderPtr->addAlignment(200536, "CAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGAC");
			alignmentReaderPtr->addAlignment(200543, "TAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAAT");
			alignmentReaderPtr->addAlignment(200552, "CCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAG");
			alignmentReaderPtr->addAlignment(200557, "GAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTC");
			alignmentReaderPtr->addAlignment(200564, "GACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGA");
			alignmentReaderPtr->addAlignment(200624, "TGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTC");
			alignmentReaderPtr->addAlignment(200647, "AAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCT");
			alignmentReaderPtr->addAlignment(200649, "ATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAA");
			alignmentReaderPtr->addAlignment(200653, "CTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAA");
			alignmentReaderPtr->addAlignment(200654, "TGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAAT");
			alignmentReaderPtr->addAlignment(200654, "TGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAAT");
			alignmentReaderPtr->addAlignment(200659, "TGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGC");
			alignmentReaderPtr->addAlignment(200668, "TCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCT");
			alignmentReaderPtr->addAlignment(200708, "ATACAATATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACGCAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACCAAGAAATTGAAGGAA");
			alignmentReaderPtr->addAlignment(200720, "GTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAAC");
			alignmentReaderPtr->addAlignment(200735, "TGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAA");
			alignmentReaderPtr->addAlignment(200750, "AAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGT");
			alignmentReaderPtr->addAlignment(200771, "CACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATAT");
			alignmentReaderPtr->addAlignment(200774, "TGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACA");
			alignmentReaderPtr->addAlignment(200775, "GAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACAC");
			alignmentReaderPtr->addAlignment(200788, "GTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATCGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAA");
			alignmentReaderPtr->addAlignment(200799, "AAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAA");
			alignmentReaderPtr->addAlignment(200805, "ATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAA");
			alignmentReaderPtr->addAlignment(200873, "CACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGG");
			alignmentReaderPtr->addAlignment(200877, "AAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGT");
			alignmentReaderPtr->addAlignment(200897, "AAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAA");
			alignmentReaderPtr->addAlignment(200903, "TAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAG");
			alignmentReaderPtr->addAlignment(200912, "CCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAGAGATCTTTT");
			alignmentReaderPtr->addAlignment(200937, "CCTAAAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTAAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAGAGATCTTTT");
			alignmentReaderPtr->addAlignment(200939, "AAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAGAGATCTTTTCCTCAGCACATGGAACATTTTGAGAGA");
			alignmentReaderPtr->addAlignment(200951, "AAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAGAGATCTTTTCCTCAGCACATGGAACATTTTGAGAGATCAGGGTAGAAA");
			alignmentReaderPtr->addAlignment(200967, "AAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAGAGAGATCTTTTCCTCAGCACATGGAACATTTTGAGAGATCAGGGTAGAAAGGTCAGAGAAATAAAA");
			alignmentReaderPtr->addAlignment(200980, "AGGGAGAGAAAAGGGTAAATAGAGAGATCTTTTCCTCAGCACATGGAACATTTTGAGAGATCAGGGTAGAAAGGTCAGAGAAATAAAAACAAACAACCTTC");
			alignmentReaderPtr->addAlignment(200980, "AGGGAGAGAAAAGGGTAAATAGAGAGATCTTTTCCTCAGCACATGGAACATTTTGAGAGATCAGGGTAGAAAGGTCAGAGAAATAAAAACAAACAACCTTC");

			referencePtr = testReferenceVariantGenerator.getReference();
			variantListPtr = testReferenceVariantGenerator.getVariants();
			alignmentReaderPtr->setRegion(referencePtr->getRegion()->getReferenceID());

		}
		*/

		static void Build31VariantsTestData(graphite::IReference::SharedPtr& referencePtr, graphite::IVariantList::SharedPtr& variantListPtr, graphite::testing::TestAlignmentReader::SharedPtr& alignmentReaderPtr)
        {
            std::string referenceSequence = "GGACCTTTCTCAGCAGCAGTGAACTTGGTAG";
            graphite::position startPosition = 0;
            std::string referenceString = std::string(startPosition, 'x');
            referenceString += referenceSequence;
			// std::cout << "asdf9" << std::endl;
            graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceString, "10", startPosition);
            testReferenceVariantGenerator.addVariant(startPosition + 1, ".", 1, {"A"});
            // testReferenceVariantGenerator.addVariant(startPosition + 3, ".", 1, {"C"});
            // testReferenceVariantGenerator.addVariant(startPosition + 20, ".", 1, {"G"});
            alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
            alignmentReaderPtr->addAlignment(2,"GACCT");
            alignmentReaderPtr->addAlignment(2,"AACCT");
            alignmentReaderPtr->addAlignment(25,"TTGGT");
            alignmentReaderPtr->addAlignment(25,"TTGGT");
            alignmentReaderPtr->addAlignment(26,"TGGTA");
            alignmentReaderPtr->addAlignment(28,"GGTAG");
            alignmentReaderPtr->addAlignment(29,"TAG");
            referencePtr = testReferenceVariantGenerator.getReference();
            variantListPtr = testReferenceVariantGenerator.getVariants();
            alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
        }

	};

	/*
	TEST(GSSWGraphManagerTests, TestBuildGraphNoVariants)
	{
		// test setup start
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		auto alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();

		graphite::position startPosition = 200000;
		std::string referenceSequence = "AAAGAGGATG";
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceSequence, "20", startPosition);
		alignmentReaderPtr->addAlignment(20000, std::string("AAAGAGGATG"));

		referencePtr = testReferenceVariantGenerator.getReference();
		variantListPtr = testReferenceVariantGenerator.getVariants();
		alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
		// test setup end

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(alignmentReaderPtr->getRegion(), 10, 5, 0);

		graphite::Variant::SharedPtr variantPtr;
		ASSERT_FALSE(variantList->getNextVariant(variantPtr));
	}

	TEST(GSSWGraphManagerTests, TestBuildGraphOneVariant)
	{
		// test setup start
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		auto alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();

		graphite::position startPosition = 200000;
		std::string referenceSequence = "AAAGAGGATG";
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceSequence, "20", startPosition);
		testReferenceVariantGenerator.addVariant(startPosition + 2, ".", 1, {"G"});
		referencePtr = testReferenceVariantGenerator.getReference();
		variantListPtr = testReferenceVariantGenerator.getVariants();
		alignmentReaderPtr->addAlignment(startPosition, std::string("AAGGAGGATG"));
		alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
		alignmentReaderManagerPtr->addAlignments(alignmentReaderPtr);
		// test setup end

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(alignmentReaderPtr->getRegion(), 20, 5, 0);

		graphite::Variant::SharedPtr variantPtr;
		ASSERT_TRUE(variantList->getNextVariant(variantPtr));

	}

	TEST(GSSWGraphManagerTests, TestBuildGraphOneVariantMismatch)
	{
		// test setup start
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		auto alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();

		graphite::position startPosition = 200000;
		std::string referenceSequence = "AAAGAGGATG";
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceSequence, "20", startPosition);
		testReferenceVariantGenerator.addVariant(startPosition + 2, ".", 1, {"G"});
		referencePtr = testReferenceVariantGenerator.getReference();
		variantListPtr = testReferenceVariantGenerator.getVariants();
		alignmentReaderPtr->addAlignment(startPosition, std::string("AAAGAG"));
		alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
		alignmentReaderManagerPtr->addAlignments(alignmentReaderPtr);
		// test setup end

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(alignmentReaderPtr->getRegion(), 20, 5, 0);

		std::map< std::string, uint32_t > alleleCount;
		alleleCount["A"] = 1;
		alleleCount["G"] = 0;


		graphite::Variant::SharedPtr variantPtr;
		ASSERT_TRUE(variantList->getNextVariant(variantPtr));
		for (auto variantAllele : variantPtr->getAlt())
		{
			ASSERT_FALSE(alleleCount.find(variantAllele) == alleleCount.end());
			ASSERT_EQ(alleleCount[variantAllele], variantPtr->getAlleleCount(variantAllele));
		}
		ASSERT_EQ(alleleCount[variantPtr->getRef()], variantPtr->getAlleleCount(variantPtr->getRef()));

	}

	TEST(GSSWGraphManagerTests, TestBuildGraphOneVariantMismatchTwoReads)
	{
		// test setup start
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		auto alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();

		graphite::position startPosition = 200000;
		std::string referenceSequence = "AAAGAGGATG";
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceSequence, "20", startPosition);
		testReferenceVariantGenerator.addVariant(startPosition + 2, ".", 1, {"G"});
		referencePtr = testReferenceVariantGenerator.getReference();
		variantListPtr = testReferenceVariantGenerator.getVariants();
		alignmentReaderPtr->addAlignment(startPosition, std::string("AAGGA"));
		alignmentReaderPtr->addAlignment(startPosition, std::string("GGATG"));
		alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
		alignmentReaderManagerPtr->addAlignments(alignmentReaderPtr);
		// test setup end

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(alignmentReaderPtr->getRegion(), 20, 5, 0);

		std::map< std::string, uint32_t > alleleCount;
		alleleCount["A"] = 0;
		alleleCount["G"] = 1;


		graphite::Variant::SharedPtr variantPtr;
		ASSERT_TRUE(variantList->getNextVariant(variantPtr));
		for (auto variantAllele : variantPtr->getAlt())
		{
			ASSERT_FALSE(alleleCount.find(variantAllele) == alleleCount.end());
			ASSERT_EQ(alleleCount[variantAllele], variantPtr->getAlleleCount(variantAllele));
		}
		ASSERT_EQ(alleleCount[variantPtr->getRef()], variantPtr->getAlleleCount(variantPtr->getRef()));

	}
*/
	/*
	TEST(GSSWGraphManagerTests, TestBuildGraphOneVariantMismatchTwoGraphs)
	{
		// test setup start
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		auto alignmentReaderPtr = std::make_shared< graphite::testing::TestAlignmentReader >();
		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();

		graphite::position startPosition = 200000;
		std::string referenceSequence = "AAAGAGGATGATGTGCAAATA";
		graphite::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceSequence, "20", startPosition);
		testReferenceVariantGenerator.addVariant(startPosition + 15, ".", 1, {"G"});
		referencePtr = testReferenceVariantGenerator.getReference();
		variantListPtr = testReferenceVariantGenerator.getVariants();
		alignmentReaderPtr->addAlignment(startPosition, std::string("AAGGA"));
		alignmentReaderPtr->addAlignment(startPosition + 10, std::string("ATGTGGAAA"));
		alignmentReaderPtr->setRegion(referencePtr->getRegion()->getRegionString());
		alignmentReaderManagerPtr->addAlignments(alignmentReaderPtr);
		// test setup end

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(alignmentReaderPtr->getRegion(), 15, 5, 2);

		std::map< std::string, uint32_t > alleleCount;
		alleleCount["C"] = 0;
		alleleCount["G"] = 1;


		graphite::Variant::SharedPtr variantPtr;
		ASSERT_TRUE(variantList->getNextVariant(variantPtr));
		for (auto variantAllele : variantPtr->getAlt())
		{
			ASSERT_FALSE(alleleCount.find(variantAllele) == alleleCount.end());
			ASSERT_EQ(alleleCount[variantAllele], variantPtr->getAlleleCount(variantAllele));
		}
		ASSERT_EQ(alleleCount[variantPtr->getRef()], variantPtr->getAlleleCount(variantPtr->getRef()));

	}
	*/

	/*
	TEST(GSSWTests, TestAlignmentReport)
	{
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		graphite::testing::TestAlignmentReader::SharedPtr alignmentReaderPtr;
		GSSWGraphTest::Build31VariantsTestData(referencePtr, variantListPtr, alignmentReaderPtr);

		auto alignmentReaderManagerPtr = std::make_shared< graphite::testing::TestAlignmentReaderManager >();
		alignmentReaderManagerPtr->addAlignments(alignmentReaderPtr);

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderManagerPtr, gsswAdjudicator);

		auto variantList = gsswGraphManager->buildGraphs(referencePtr->getRegion(), 3000, 1000, 100);

		// std::ofstream outStream;
		// outStream.open("og.txt", std::ios::out);
		// graphite::gssw::AlignmentReporter::Instance()->printAlignmentReportsToStream(outStream);
		// outStream.close();

		std::ofstream outVCF;
		outVCF.open("test.vcf", std::ios::out);
		variantList->printToVCF(outVCF);
		outVCF.close();
	}
	*/

	/*
	TEST(GSSWTests, TestConstructChr20)
	{
		std::string fastaPath = TEST_FASTA_FILE;
		// std::string vcfPath = TEST_SMALL_MEI_VCF_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;
		std::string bamPath = TEST_BAM_FILE;

		graphite::Region::SharedPtr regionPtr = std::make_shared< graphite::Region >("20");
		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReaderPtr = std::make_shared< graphite::VariantListVCFPreloaded >(vcfPath, regionPtr);
		// auto bamAlignmentReaderManager = std::make_shared< graphite::BamAlignmentReaderManager >(bamPath);
		// auto graphManagerPtr = std::make_shared< graphite::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderManager, 25);
		auto bamAlignmentReaderPreloadManager = std::make_shared< graphite::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);
		vcfFileReaderPtr->loadVariantsFromFile();

		std::cout << "Finished loading BAM and VCF" << std::endl;

		auto gsswAdjudicator = std::make_shared< graphite::gssw::GSSWAdjudicator >();
		auto gsswGraphManager = std::make_shared< graphite::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderPreloadManager, gsswAdjudicator);
		auto variantList = gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

		std::cout << "starting to print vcf" << std::endl;

		std::ofstream outVCF;
		outVCF.open("output.vcf", std::ios::out);
		variantList->printToVCF(outVCF);
		outVCF.close();

		std::cout << "printing other output" << std::endl;

		std::ofstream outStream;
		outStream.open("og.txt", std::ios::out);
		graphite::gssw::AlignmentReporter::Instance()->printAlignmentReportsToStream(outStream);
		outStream.close();
	}
	*/

	TEST(GSSWTests, TestConstructTestData)
	{
		/*
		graphite::IReference::SharedPtr referencePtr;
		graphite::IVariantList::SharedPtr variantListPtr;
		graphite::testing::TestAlignmentReader::SharedPtr alignmentReaderPtr;
		GSSWGraphTest::Build31VariantsTestData(referencePtr, variantListPtr, alignmentReaderPtr);

		auto gsswGraph = std::make_shared< GSSWGraphTest >(referencePtr, variantListPtr, alignmentReaderPtr);
		gsswGraph->constructGraph();
		*/
	}



}

#endif //GRAPHITE_GSSW_GSSWGRAPH_TESTS_HPP
