#ifndef GRAPHITE_GRAPHPATHREFTESTS_HPP
#define GRAPHITE_GRAPHPATHREFTESTS_HPP

#include "TestConfig.h"
#include "core/region/Region.h"
#include "core/reference/Reference.h"
#include "core/reference/FastaReference.h"
#include "core/graph/GSSWGraph.h"

#include "core/file/FastaFileWriter.h"

TEST(GraphPathRefTests, ConstructPaths)
{
    uint32_t readLength = 5;  
    std::string vcfLine = "1\t10\trs11575e97\tT\tC\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line.

	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);

    auto variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), referencePtr, readLength);
    std::vector< graphite::IVariant::SharedPtr > variantPtrs = {variantPtr};
    auto variantListPtr = std::make_shared< graphite::VariantList >(variantPtrs, referencePtr);
    
    auto gsswRegionPtr = std::make_shared< graphite::Region >("1", variantPtr->getPosition() - readLength, variantPtr->getPosition() + readLength, graphite::Region::BASED::ONE);
    //auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, regionPtr, 1, 1, 1, 1, numGraphCopies);
    auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, gsswRegionPtr, 1, 1, 1, 1, 1);
    gsswGraphPtr->constructGraph();
    
    // Need to set up a graphManager to extract the graph path headers and sequences in the above code. After that, replacing the gsswGraphPtr with a graphManager object should fix the rest of the code.
    /*
    std::vector< std::string > graphPathHeaders = gsswGraphPtr->getGraphPathHeaders();
    std::vector< std::string > graphPathSequences = gsswGraphPtr->getGraphPathSequences();
    std::vector< int > graphPathLengths = gsswGraphPtr->getGraphPathLengths();

    std::cout << "Number of headers: " << graphPathHeaders.size() << std::endl;
    std::cout << "Number of sequences: " << graphPathSequences.size() << std::endl;
    for (int i = 0; i < graphPathHeaders.size() && i < graphPathSequences.size(); ++i)
    {
        std::cout << graphPathHeaders[i] << "\t" << "Length: " << graphPathLengths[i] << std::endl;
        std::cout << graphPathSequences[i] <<  std::endl;
    }
    */
}

// May not need this test since the previous test already prints the headers and sequences to screen.
/*
TEST(GraphPathRefTests, WritePathsToFasta)
{
    uint32_t readLength = 5;  
    std::string vcfLine = "1\t10\trs11575e97\tT\tC\t34439.5\tPASS\tAA=G;AC=22;AF=0.0178427;AN=1233;DP=84761;NS=1233;AMR_AF=0.0000;AFR_AF=0.0000;EUR_AF=0.0000;SAS_AF=0.0000;EAS_AF=0.0451\tGT\t0\t0"; // is not the complete first line.

	auto fastaRegionPtr = std::make_shared< graphite::Region >("1", graphite::Region::BASED::ONE);
	auto referencePtr = std::make_shared< graphite::FastaReference >(TEST_FASTA_FILE, fastaRegionPtr);

    auto variantPtr = graphite::Variant::BuildVariant(vcfLine.c_str(), referencePtr, readLength);
    std::vector< graphite::IVariant::SharedPtr > variantPtrs = {variantPtr};
    auto variantListPtr = std::make_shared< graphite::VariantList >(variantPtrs, referencePtr);
    
    auto gsswRegionPtr = std::make_shared< graphite::Region >("1", variantPtr->getPosition() - readLength, variantPtr->getPosition() + readLength, graphite::Region::BASED::ONE);
    //auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, regionPtr, 1, 1, 1, 1, numGraphCopies);
    auto gsswGraphPtr = std::make_shared< graphite::GSSWGraph >(referencePtr, variantListPtr, gsswRegionPtr, 1, 1, 1, 1, 1);
    gsswGraphPtr->constructGraph();
    
    std::vector< std::string > graphPathHeaders = gsswGraphPtr->getGraphPathHeaders();
    std::vector< std::string > graphPathSequences = gsswGraphPtr->getGraphPathSequences();

    graphite::FastaFileWriter fastaFileWriter;
    fastaFileWriter.open("");
    fastaFileWriter.write(graphPathHeaders, graphPathSequences);
    fastaFileWriter.close();
}
*/

#endif // GRAPHITE_GRAPHPATHREFTESTS_HPP
