#ifndef GWIZ_GSSW_GSSWGRAPH_TESTS_HPP
#define GWIZ_GSSW_GSSWGRAPH_TESTS_HPP

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
#include "plugins/gssw/graph/GSSWGraph.h"
#include "plugins/gssw/graph/AlignmentReporter.h"

#include "core/alignments/BamAlignmentReader.h"
#include "core/alignments/BamAlignmentReaderManager.h"
#include "core/alignments/BamAlignmentReaderPreloadManager.h"
#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/reference/FastaReference.h"
#include "core/utils/ThreadPool.hpp"


namespace
{

	class GSSWGraphTest : public gwiz::gssw::GSSWGraph
	{
	public:
		typedef std::shared_ptr< GSSWGraphTest > GSSWGraphTestPtr;

		GSSWGraphTest(gwiz::IReference::SharedPtr referencePtr, gwiz::IVariantList::SharedPtr variantListPtr, gwiz::IAlignmentReader::SharedPtr alignmentReader) :
			gwiz::gssw::GSSWGraph(referencePtr, variantListPtr, alignmentReader)
		{
		}

		~GSSWGraphTest()
		{
		}

		/*
		gwiz::gssw::GSSWGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(gwiz::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		gwiz::gssw::GSSWGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}
		*/

		static void Build31VariantsTestData(gwiz::IReference::SharedPtr& referencePtr, gwiz::IVariantList::SharedPtr& variantListPtr, gwiz::testing::TestAlignmentReader::SharedPtr& alignmentReaderPtr)
		{
			std::string referenceSequence = "ATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAG";

			gwiz::position startPosition = 200000;
			std::string referenceString = std::string(startPosition, 'x');
			referenceString += referenceSequence;
			gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceString, "20", startPosition);
            testReferenceVariantGenerator.addVariant(startPosition + 14, ".", 1, {"G"});
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

			alignmentReaderPtr = std::make_shared< gwiz::testing::TestAlignmentReader >();
			alignmentReaderPtr->addAlignment(199924, "GGACCTTTCTCAGCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCCCCAACTCTAGG");
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
			alignmentReaderPtr->setRegion(referencePtr->getRegion());

			/*
			std::cout << "start position: " << referencePtr->getRegion()->getStartPosition() << std::endl;
			std::cout << "end position: " << referencePtr->getRegion()->getEndPosition() << std::endl;
			std::cout << "ref len: " << referencePtr->getRegion()->getEndPosition() - referencePtr->getRegion()->getStartPosition() << std::endl;
			std::cout << "length: " << referenceString.size() << std::endl;
			*/
		}

	};

	TEST(GSSWTests, TestAlignmentReport)
	{
		gwiz::IReference::SharedPtr referencePtr;
		gwiz::IVariantList::SharedPtr variantListPtr;
		gwiz::testing::TestAlignmentReader::SharedPtr alignmentReaderPtr;
		GSSWGraphTest::Build31VariantsTestData(referencePtr, variantListPtr, alignmentReaderPtr);

		auto gsswGraph = std::make_shared< gwiz::gssw::GSSWGraph >(referencePtr, variantListPtr, alignmentReaderPtr);
		gsswGraph->constructGraph();

		gwiz::gssw::AlignmentReporter::Instance()->printAlignmentReportsToStream(std::cout);
		// auto graphManagerPtr = std::make_shared< gwiz::gssw::GraphManager >(referencePtr, variantListPtr, alignmentReaderPtr, 25);
		// graphManagerPtr->buildGraphs();
		// auto gsswGraph = std::make_shared< GSSWGraphTest >(referencePtr, variantListPtr, alignmentReaderPtr);
		// gsswGraph->constructGraph();
	}

	TEST(GSSWTests, TestConstructChr20)
	{
		// boost::function< void, (int) > funct [](int x) { std::cout << "pool: " << x << std::endl; };
		/*
		std::function<void (int*)> f2 = [](int* x){ *x += 1; std::this_thread::sleep_for(std::chrono::seconds(1)); };

		int count = 0;
		int count2 = 0;
		for (int i = 0; i < 100; ++i)
		{
			++count2;
			gwiz::ThreadPool::Instance()->postJob(boost::bind(f2, &count));
		}
		std::this_thread::sleep_for(std::chrono::seconds(60));
		// gwiz::ThreadPool::Instance()->joinAll();
		std::cout << "Count: " << count << std::endl;
		std::cout << "Count2: " << count2 << std::endl;

		if (true) { return; }
		*/
		// IReference::SharedPtr referencePtr, IVariantList::SharedPtr variantListPtr, IAlignmentReaderManager::SharedPtr alignmentReaderManager, size_t padding) :

		/*
		std::string fastaPath = TEST_FASTA_FILE;
		std::string vcfPath = TEST_1KG_CHR20_VCF_FILE;
		std::string bamPath = TEST_BAM_FILE;
		// bamPath = "/d1/data/project_bam/platinum_genomes/NA12878_2_sorted_Dedup_realign_recal.bam";
		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20");
		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReaderPtr = std::make_shared<gwiz::VCFFileReader>(vcfPath);
		// auto bamAlignmentReaderManager = std::make_shared< gwiz::BamAlignmentReaderManager >(bamPath);
		// auto graphManagerPtr = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderManager, 25);
		auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);
		auto graphManagerPtr = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderPreloadManager, 25);
		graphManagerPtr->buildGraphs();
		*/
	}

	TEST(GSSWTests, TestConstructTestData)
	{
		/*
		gwiz::IReference::SharedPtr referencePtr;
		gwiz::IVariantList::SharedPtr variantListPtr;
		gwiz::testing::TestAlignmentReader::SharedPtr alignmentReaderPtr;
		GSSWGraphTest::Build31VariantsTestData(referencePtr, variantListPtr, alignmentReaderPtr);

		auto gsswGraph = std::make_shared< GSSWGraphTest >(referencePtr, variantListPtr, alignmentReaderPtr);
		gsswGraph->constructGraph();
		*/

/*
		std::cout << "contigs constructed" << std::endl;

		auto contigs = variantGraph->getContigs();
		size_t contigCount = 0;
		while (!contigs.empty())
		{
			auto contig = contigs.front();
			contigCount += contig->getContigs().size();
			contigs.pop();
		}
		std::cout << "Contig Count: " << contigCount << std::endl;
*/

		// variantGraph->printGraph("test1.dot");


		// gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(m_reference_string, "20", 6000);
		// testReferenceVariantGenerator.addVariant(6010, ".", {"A","C"});
		// auto variantGraph = std::make_shared< gwiz::vg::VariantGraph >(testReferenceVariantGenerator.getReference(), testReferenceVariantGenerator.getVariants());

	}

}

#endif //GWIZ_GSSW_GSSWGRAPH_TESTS_HPP
