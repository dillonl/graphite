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


		gwiz::gssw::GSSWGraph::VariantVertexDescriptor getReferenceVertexContainsPositionTest(gwiz::position pos)
		{
			return getReferenceVertexContainsPosition(pos);
		}

		gwiz::gssw::GSSWGraph::GraphPtr getGraph()
		{
			return this->m_graph_ptr;
		}


		static void Build31VariantsTestData(gwiz::IReference::SharedPtr& referencePtr, gwiz::IVariantList::SharedPtr& variantListPtr, gwiz::testing::TestAlignmentReader::SharedPtr& alignmentReaderPtr)
		{
			std::string referenceSequence = "ATGTCACCTCTCCCCCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGATGTAGAGGGAGAGAAAAGGGTAAATAG";

			gwiz::position startPosition = 200000;
			std::string referenceString = std::string(startPosition, 'x');
			referenceString += referenceSequence;
			gwiz::testing::TestReferenceVariantGenerator testReferenceVariantGenerator(referenceString, "20", startPosition);
            testReferenceVariantGenerator.addVariant(startPosition + 14, ".", 1, {"GTTTACTATCGTACAGTACGATACAGTACAGTACCCGTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAA"});
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
			alignmentReaderPtr->addAlignment(199924, "GGACCTTTCTCAGCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGT");
			alignmentReaderPtr->addAlignment(199933, "TCAGCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGA");
			alignmentReaderPtr->addAlignment(199936, "GCAGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGATAC");
			alignmentReaderPtr->addAlignment(199938, "AGCAGTGAACTTGGGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGATACAG");
			alignmentReaderPtr->addAlignment(199951, "GGGTGCTCACAACCTGTGCAAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGATACAGTACAGTACCCGTT");
			alignmentReaderPtr->addAlignment(199970, "AAAACCAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGATACAGTACAGTACCCGTTGGTGACACACGTCGACATG");
			alignmentReaderPtr->addAlignment(199975, "CAGCTGTGGTGGCTAAGGATTGCCTATGTCACCTCTCCGTTTACTATCGTACAGTACGATACAGTACAGTACCCGTTGGTGACACACGTCGACATGCTACA");
			alignmentReaderPtr->addAlignment(200004, "CACCTCTCCGTTTACTATCGTACAGTACGATACAGTACAGTACCCGTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGC");
			alignmentReaderPtr->addAlignment(200017, "ACTATCGTACAGTACGATACAGTACAGTACCCGTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATA");
			alignmentReaderPtr->addAlignment(200042, "TACCCGTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGC");
			alignmentReaderPtr->addAlignment(200046, "CGTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATG");
			alignmentReaderPtr->addAlignment(200047, "GTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGC");
			alignmentReaderPtr->addAlignment(200047, "GTTGGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGC");
			alignmentReaderPtr->addAlignment(200050, "GGTGACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGCAGC");
			alignmentReaderPtr->addAlignment(200054, "ACACACGTCGACATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGCAGCTTGG");
			alignmentReaderPtr->addAlignment(200066, "ATGCTACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCT");
			alignmentReaderPtr->addAlignment(200070, "TACAGCTACATCGACCGTAGCCACACAGCTCGCGTGTGTGTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCA");
			alignmentReaderPtr->addAlignment(200109, "GTATATATATAAAAAAAAAAAACCAACTCTAGGCAATGCAGCTTGGGGATAGACTCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCC");
			alignmentReaderPtr->addAlignment(200163, "TCCTTCCACTTGGGGGAAGAAGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCC");
			alignmentReaderPtr->addAlignment(200183, "AGAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGT");
			alignmentReaderPtr->addAlignment(200184, "GAGGGAAGAGTACAGAGGGCTTTGCCTTGCAACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTG");
			alignmentReaderPtr->addAlignment(200214, "AACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCA");
			alignmentReaderPtr->addAlignment(200215, "ACTTGGGTACCAGCTCAGCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAA");
			alignmentReaderPtr->addAlignment(200232, "GCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGA");
			alignmentReaderPtr->addAlignment(200232, "GCCACAGTAAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGA");
			alignmentReaderPtr->addAlignment(200240, "AAAGTAAAGTATCAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACT");
			alignmentReaderPtr->addAlignment(200252, "CAAAAGTTACCCAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGG");
			alignmentReaderPtr->addAlignment(200263, "CAGCATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAG");
			alignmentReaderPtr->addAlignment(200267, "ATGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAA");
			alignmentReaderPtr->addAlignment(200268, "TGGTGCCAGCTGTGGTGGCCACTGGACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAG");
			alignmentReaderPtr->addAlignment(200292, "GACTTGCCCTTCCCCCAACTCCAAGCAGCCTGGCACAGAGAGAGAGACTCCTTTTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAA");
			alignmentReaderPtr->addAlignment(200345, "TTGTTTGGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAG");
			alignmentReaderPtr->addAlignment(200351, "GGGGGTAAATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCA");
			alignmentReaderPtr->addAlignment(200358, "AATGAGGGAAGAGAAGAAGAAACTCTGCCTGGTAACCCAGGGAATTTGGCCAAATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAG");
			alignmentReaderPtr->addAlignment(200410, "AATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATA");
			alignmentReaderPtr->addAlignment(200410, "AATTTAAACCCCAGCCCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATA");
			alignmentReaderPtr->addAlignment(200425, "CCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAAT");
			alignmentReaderPtr->addAlignment(200425, "CCACTAAGGTGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAAT");
			alignmentReaderPtr->addAlignment(200434, "TGGTTCCTCTAGGACTCAGCGAGAGTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACT");
			alignmentReaderPtr->addAlignment(200458, "GTTGCAGTGTTTCTGAGCTTAGGGCACCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAA");
			alignmentReaderPtr->addAlignment(200484, "CCCTCTAGTGCTGATATAGTTTCAATAATCACAGGCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAA");
			alignmentReaderPtr->addAlignment(200517, "GCTCAAATCACAACACTCAATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCT");
			alignmentReaderPtr->addAlignment(200536, "ATCTCCTTCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGAC");
			alignmentReaderPtr->addAlignment(200543, "TCAAATACCTGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACA");
			alignmentReaderPtr->addAlignment(200552, "TGAAAAGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTC");
			alignmentReaderPtr->addAlignment(200557, "AGCCTTCCCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAG");
			alignmentReaderPtr->addAlignment(200564, "CCAAGAAGGATGGGTGCAAACAAGCCCAGATTGTGAAGGCTACAATATGTATCTAACTCTTCAATGCCCAGACATCAACAACCATCTTCAAGAGTTAAGAA");
			alignmentReaderPtr->addAlignment(200624, "GACATCAACAACCATCTTCAAGAGTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAG");
			alignmentReaderPtr->addAlignment(200647, "GTTAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATA");
			alignmentReaderPtr->addAlignment(200649, "TAAGAACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAAT");
			alignmentReaderPtr->addAlignment(200653, "AACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAA");
			alignmentReaderPtr->addAlignment(200654, "ACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAA");
			alignmentReaderPtr->addAlignment(200654, "ACATCCAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAA");
			alignmentReaderPtr->addAlignment(200659, "CAGGGAAATATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGC");
			alignmentReaderPtr->addAlignment(200668, "ATGACCTCATCAAATGAACTAAATAAGGCATCAGTGACCAATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGAT");
			alignmentReaderPtr->addAlignment(200708, "ATCTGAGAATGATGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTC");
			alignmentReaderPtr->addAlignment(200720, "TGGAGATATGTGACTTTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCA");
			alignmentReaderPtr->addAlignment(200735, "TTTTAGACAAATAATTCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAA");
			alignmentReaderPtr->addAlignment(200750, "TCAAAATAGCTGTCTTGATGAAGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATT");
			alignmentReaderPtr->addAlignment(200771, "AGCTCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAAT");
			alignmentReaderPtr->addAlignment(200774, "TCAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCT");
			alignmentReaderPtr->addAlignment(200775, "CAACAAACTTCAAGACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTG");
			alignmentReaderPtr->addAlignment(200788, "GACAACACAGAGAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTT");
			alignmentReaderPtr->addAlignment(200799, "GAAAGAATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAA");
			alignmentReaderPtr->addAlignment(200805, "ATTCAGAATTACATCAGAGAAGTTTCACAAAGAAATTGAAGTAATTTCTAAAAAATCAAGCAGAAATTCTGGCACTGAAAAGTTTGATTGTCAAAGTGAAA");
			alignmentReaderPtr->addAlignment(200873, "CTGGCACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAA");
			alignmentReaderPtr->addAlignment(200877, "CACTGAAAAGTTTGATTGTCAAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATAT");
			alignmentReaderPtr->addAlignment(200897, "AAAGTGAAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAA");
			alignmentReaderPtr->addAlignment(200903, "AAAAATGCATAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAA");
			alignmentReaderPtr->addAlignment(200912, "TAAGAGTCTTTCAACGGCAGAATTGATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAA");
			alignmentReaderPtr->addAlignment(200937, "ATCAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATA");
			alignmentReaderPtr->addAlignment(200939, "CAAGCAGAAGGAACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGT");
			alignmentReaderPtr->addAlignment(200951, "ACTGGTGAGAACTGGCTATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTA");
			alignmentReaderPtr->addAlignment(200967, "TATCCAAATATACACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTG");
			alignmentReaderPtr->addAlignment(200980, "ACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGA");
			alignmentReaderPtr->addAlignment(200980, "ACGAAGCCAAAAAAAGAATTAAAAAAGAATAAAGTATGCCTACAAAATGTAGAAAATAGTCTCAAAAGGGTAAATCCAAGAGTTATTGGTCTTAAAGAGGA");

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
/*
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
*/
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
		std::ofstream ofs ("test.txt",std::ofstream::out);
		ofs << "hello";
		// bamPath = "/d1/data/project_bam/platinum_genomes/NA12878_2_sorted_Dedup_realign_recal.bam";
/*
		gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >("20");
		auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
		auto vcfFileReaderPtr = std::make_shared<gwiz::VCFFileReader>(vcfPath);
		// auto bamAlignmentReaderManager = std::make_shared< gwiz::BamAlignmentReaderManager >(bamPath);
		// auto graphManagerPtr = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderManager, 25);
		auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);
		auto graphManagerPtr = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderPreloadManager, 25);
		graphManagerPtr->buildGraphs();
		gwiz::gssw::AlignmentReporter::Instance()->printAlignmentReportsToStream(ofs);

		ofs.close();
		std::cout << "We are here" << std::endl;
	}
*/
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
