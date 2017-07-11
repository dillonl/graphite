#include "core/alignment/AlignmentManager.hpp"
#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/variant/VCFFileReader.h"
#include "core/reference/FastaReference.h"
#include "core/mapping/MappingManager.h"
#include "core/variant/VCFHeader.h"
#include "core/util/Params.h"
#include "core/util/ThreadPool.hpp"
#include "core/graph/GraphManager.h"
#include "core/adjudicator/GSSWAdjudicator.h"
#include "core/variant/VCFHeader.h"
#include "core/sample/SampleManager.h"
#include "core/region/Region.h"
#include "core/file/IFile.h"
#include "core/file/BGZFFileWriter.h"
#include "core/file/ASCIIFileWriter.h"
#include "core/file/BamHeaderReader.h"

#include "core/file/FastaFileWriter.h"
#include "core/util/Utility.h"

#include <thread>
#include <unordered_set>
#include <fstream>
#include <stdio.h>

#include <unordered_map>
#include <string>

#include <zlib.h>   // May not need this.

int main(int argc, char** argv)
{
	unsigned long milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	graphite::Params params;
	params.parseGSSW(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		exit(0);
	}
	auto bamPaths = params.getBAMPaths();
	auto fastaPath = params.getFastaPath();
	auto vcfPaths = params.getInVCFPaths();
	auto outputDirectory = params.getOutputDirectory();
	auto paramRegionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto excludeDuplicates = params.getExcludeDuplicates();
	auto graphSize = params.getGraphSize();
	graphite::FileType fileType = graphite::FileType::ASCII;

	graphite::ThreadPool::Instance()->setThreadCount(threadCount);

	std::vector< graphite::Region::SharedPtr > regionPtrs;
	if (paramRegionPtr == nullptr)
	{
		regionPtrs = graphite::VCFFileReader::GetAllRegionsInVCF(vcfPaths);
	}
	else
	{
		regionPtrs.emplace_back(paramRegionPtr);
	}

	uint32_t readLength = graphite::BamAlignmentManager::GetReadLength(bamPaths);

    // One solution for multiple bam files is to make a struct or class that contains all the open new_bam files. Then use them to write out the BamAlignment to the appropriate new_bam file.
    // Retrieve SAM header from input bam file.
    /*
    if (bamPaths.size() == 1)
    {
        graphite::BamAlignmentReader bar(bamPaths[0]);
        bar.open();
        std::string samHeader = bar.getSamHeader();
        bar.close();
        std::ofstream samFile;
        samFile.open("NewSamFile.sam", std::ios::trunc);
        samFile << samHeader;
        samFile.close();
    }
    */
    /*
    for (auto& bamPath : bamPaths)
    {
        graphite::BamAlignmentReader bar(bamPath);
        std::string samHeader = bar.getSamHeader();
        std::fstream samFile;
        samFile.open("NewSamFile.sam");//, std::ios::trunc);
        samFile << samHeader << std::endl;
        samFile.close();
    }
    */

	graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(graphite::BamAlignmentManager::GetSamplePtrs(bamPaths));

	auto alignmentReaderManagerPtr = std::make_shared< graphite::AlignmentReaderManager< graphite::BamAlignmentReader > >(bamPaths, threadCount);

	std::unordered_map< std::string, graphite::IFileWriter::SharedPtr > vcfoutPaths;
	for (auto vcfPath : vcfPaths)
	{
		std::string path = vcfPath.substr(vcfPath.find_last_of("/") + 1);
		std::string filePath = outputDirectory + "/" + path;
		uint32_t counter = 1;
		while (graphite::IFile::fileExists(filePath, false))
		{
			std::string extension = vcfPath.substr(vcfPath.find_last_of(".") + 1);
			std::string fileNameWithoutExtension = path.substr(0, path.find_last_of("."));
			filePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + extension;
			++counter;
		}
		// filePath += ".tmp";
		graphite::IFileWriter::SharedPtr fileWriterPtr;
		if (fileType == graphite::FileType::BGZF)
		{
			fileWriterPtr = std::make_shared< graphite::BGZFFileWriter >(filePath);
		}
		else
		{
			fileWriterPtr = std::make_shared< graphite::ASCIIFileWriter >(filePath);
		}
		fileWriterPtr->open();
		vcfoutPaths[vcfPath] = fileWriterPtr;
	}

	std::unordered_set< std::string > outputPaths;
	bool firstTime = true;

    // Refactor code so that all information related to the graphPathHeaders and the resulting bam can be written out after this for loop.
    // GraphPaths fasta file map
    std::unordered_map< std::string, std::string > headerSequenceMap; 

    // GraphPaths bed file map
    std::unordered_map< uint32_t, graphite::NodeInfo::SharedPtr > nodeInfoMap;

    // Create a headerSequenceMap.
    // Add code inside for loop that will insert key:values into the map.
    // Test to see if I'm appropriately writing out the graphPaths and seuqences.
    // If it works then I can mark the graphPathHeader map in the GraphManager to be removed. Need to implement the BED and SAM stuff before actually deleting it.

	for (uint32_t regionCount = 0; regionCount < regionPtrs.size(); ++regionCount)
	{
		auto regionPtr = regionPtrs[regionCount];

		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

		// load variants from vcf
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, readLength);
		variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously

		variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory

		// load bam alignments
		auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(sampleManagerPtr, regionPtr, alignmentReaderManagerPtr, excludeDuplicates);
		bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
		bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory


		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
		bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

		std::deque< std::shared_ptr< std::future< void > > > variantManagerFutureFunctions;
		for (auto& iter : variantManagerPtr->getVCFReadersAndVariantListsMap())
		{
			auto futureFunct = graphite::ThreadPool::Instance()->enqueue(std::bind(&graphite::IVariantList::processOverlappingAlleles, iter.second));
			variantManagerFutureFunctions.push_back(futureFunct);
		}
		while (!variantManagerFutureFunctions.empty())
		{
			variantManagerFutureFunctions.front()->wait();
			variantManagerFutureFunctions.pop_front();
		}

		milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

		// create an adjudicator for the graph
		auto gsswAdjudicator = std::make_shared< graphite::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);

		// the gsswGraphManager adjudicates on the variantManager's variants
		auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
		// auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, alignmentManager, gsswAdjudicator);
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength);

        // Append headers and sequences to headerSequenceMap
        std::unordered_map< std::string, std::string > tempHeaderSequenceMap = gsswGraphManager->getHeaderSequenceMap();
        headerSequenceMap.insert(tempHeaderSequenceMap.begin(), tempHeaderSequenceMap.end());

        // Store BED entries for BED file.
        std::unordered_map< uint32_t, graphite::NodeInfo::SharedPtr > tempNodeInfoMap = gsswGraphManager->getNodeInfoMap();
        nodeInfoMap.insert(tempNodeInfoMap.begin(), tempNodeInfoMap.end());

		graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);
		graphite::MappingManager::Instance()->clearRegisteredMappings();

		std::vector< std::shared_ptr< std::thread > > fileWriters;
		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFReadersAndVariantListsMap();
		std::deque< std::shared_ptr< std::future< void > > > vcfWriterFutureFunctions;
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			auto vcfReaderPtr = iter.first;
			auto vcfPath = vcfReaderPtr->getFilePath();
			graphite::IFileWriter::SharedPtr fileWriter = vcfoutPaths[vcfPath];
			std::string currentVCFOutPath = fileWriter->getFilePath();
			auto variantListPtr = iter.second;
			auto vcfHeaderPtr = vcfReaderPtr->getVCFHeader();

			vcfHeaderPtr->registerActiveSample(sampleManagerPtr);
			if (firstTime)
			{
				outputPaths.emplace(currentVCFOutPath);
			}
			auto funct = std::bind(&graphite::VariantList::writeVariantList, variantListPtr, fileWriter, vcfHeaderPtr, firstTime);
			auto functFuture = graphite::ThreadPool::Instance()->enqueue(funct);
			vcfWriterFutureFunctions.push_back(functFuture);
		}

		while (!vcfWriterFutureFunctions.empty())
		{
			vcfWriterFutureFunctions.front()->wait();
			vcfWriterFutureFunctions.pop_front();
		}

		firstTime = false;
	}

    // Write out GraphPaths fasta.
    graphite::FastaFileWriter fastaWriter;
    fastaWriter.open("GraphPaths.fa");
    for (auto& hs: headerSequenceMap)
    {
        fastaWriter.write(hs.first, hs.second);
    }
    fastaWriter.close();

    // Write out GraphPaths bed.
    std::ofstream bedFile;
    bedFile.open("GraphPaths.bed");
    for (auto& hs: headerSequenceMap)
    {
        uint8_t nodeStringStart= hs.first.find("_") + 1;
        std::string nodeString = hs.first.substr(nodeStringStart);
        std::vector< std::string > nodeVector;
        graphite::split(nodeString, ':', nodeVector);

        int32_t startPosition = 0;
        int32_t endPosition;
        std::string variantType;
        for (int i = 0; i < nodeVector.size(); ++i)
        {
            auto iter = nodeInfoMap.find(std::stoi(nodeVector[i]));
            //std::unordered_map< uint32_t, graphite::NodeInfo::SharedPtr >::const_iterator iter = nodeInfoMap.find("2");
            endPosition = startPosition + iter->second->getLength();
            if (iter->second->getVariantType() == 0)
                variantType = "Ref";
            else
                variantType = "Alt";

            bedFile 
                << hs.first << "\t" 
                << startPosition << "\t" 
                << endPosition << "\t"
                << variantType
                << std::endl;
            startPosition = endPosition;
        }
    }
    bedFile.close();
    
    // Write out SAM header and alignments.
    graphite::BamHeaderReader bamFile(bamPaths[0]);
    bamFile.open();
    for (auto& hs: headerSequenceMap)
    {
        bamFile.addPathHeaderToSamHeader(hs.first, hs.second.length());
    }
    bamFile.addReadGroupsToSamHeader();
    std::string samHeader = bamFile.getModifiedSamHeader();
    bamFile.close();

    std::ifstream tempAlignmentFile("TempAlignmentFile.sam");
    std::ofstream samFile("NewSamFile.sam", std::ios::app);
    samFile << samHeader;
    for (std::string str; std::getline(tempAlignmentFile, str); )
    {
        samFile << str << std::endl;
    }
    samFile.close();
    tempAlignmentFile.close();

    // Remove temporary file.
    remove("TempAlignmentFile.sam");
    
    // Write out vcfs
	for (auto& iter : vcfoutPaths)
	{
		graphite::IFileWriter::SharedPtr fileWriter = iter.second;
		fileWriter->close();
	}

	// graphite::GSSWAdjudicator* adj_p;
	// std::cout << "adj counts: " << (uint32_t)adj_p->s_adj_count << " [total]" << std::endl;

	return 0;
}
