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
#include "core/file/SAMFileWriter.h"

#include "core/util/Utility.h"

#include <thread>
#include <unordered_set>
#include <fstream>
#include <stdio.h>

#include <unordered_map>
#include <string>


void updateFileMap (std::unordered_map< std::string, graphite::IFileWriter::SharedPtr > &outputFileMap, std::string inputFilePath, graphite::FileType fileType, std::string outputDirectory, std::string fileExtension)
{
    std::string fileName = inputFilePath.substr(inputFilePath.find_last_of("/") + 1);
    std::string fileNameWithoutExtension = fileName.substr(0, fileName.find_last_of("."));
    std::string outputFilePath = outputDirectory + "/" + fileNameWithoutExtension + "." +  fileExtension;
    uint32_t counter = 1;
    while (graphite::IFile::fileExists(outputFilePath, false))
    {
        outputFilePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + fileExtension;
        ++counter;
    }
    graphite::IFileWriter::SharedPtr fileWriterPtr;
    if (fileType == graphite::FileType::BGZF)
        fileWriterPtr = std::make_shared< graphite::BGZFFileWriter >(outputFilePath);
    else if (fileType == graphite::FileType::ASCII)
        fileWriterPtr = std::make_shared< graphite::ASCIIFileWriter >(outputFilePath);
    else
        fileWriterPtr = std::make_shared< graphite::SAMFileWriter >(outputFilePath);

    fileWriterPtr->open();
    outputFileMap[inputFilePath] = fileWriterPtr;
}

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
    auto isIGVOutput = params.getIGVOutput();

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

	graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(graphite::BamAlignmentManager::GetSamplePtrs(bamPaths));

	auto alignmentReaderManagerPtr = std::make_shared< graphite::AlignmentReaderManager< graphite::BamAlignmentReader > >(bamPaths, threadCount);

    // Create file map for tracking file writers.
	std::unordered_map< std::string, graphite::IFileWriter::SharedPtr > outputFileMap;
    graphite::FileType asciiFileType = graphite::FileType::ASCII;
    for (auto& vcfPath : vcfPaths)
    {
        updateFileMap (outputFileMap, vcfPath, asciiFileType, outputDirectory, "vcf");
    }

    // Append temporary SAM file to outputFileMap.
    std::string firstFileName_withoutExtension = vcfPaths[0].substr(0, vcfPaths[0].find_last_of("."));
    std::string tempSamFileName = firstFileName_withoutExtension + "TEMP" + "." + "sam";
    graphite::FileType samFileType = graphite::FileType::SAM;
    updateFileMap (outputFileMap, tempSamFileName, samFileType, outputDirectory, "sam");
    graphite::SAMFileWriter::SharedPtr tempSamFilePtr = std::dynamic_pointer_cast< graphite::SAMFileWriter >(outputFileMap.at(tempSamFileName));

	std::unordered_set< std::string > outputPaths;
	bool firstTime = true;

    // GraphPaths fasta file map
    std::unordered_map< std::string, std::string > headerSequenceMap; 

    // GraphPaths bed file map
    std::unordered_map< uint32_t, graphite::NodeInfo::SharedPtr > nodeInfoMap;

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
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength, isIGVOutput, tempSamFilePtr);

        if (isIGVOutput)
        {
            // Append headers and sequences to headerSequenceMap
            std::unordered_map< std::string, std::string > tempHeaderSequenceMap = gsswGraphManager->getHeaderSequenceMap();
            headerSequenceMap.insert(tempHeaderSequenceMap.begin(), tempHeaderSequenceMap.end());

            // Store BED entries for BED file.
            std::unordered_map< uint32_t, graphite::NodeInfo::SharedPtr > tempNodeInfoMap = gsswGraphManager->getNodeInfoMap();
            nodeInfoMap.insert(tempNodeInfoMap.begin(), tempNodeInfoMap.end());
        }

		graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);
		graphite::MappingManager::Instance()->clearRegisteredMappings();

		std::vector< std::shared_ptr< std::thread > > fileWriters;
		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFReadersAndVariantListsMap();
		std::deque< std::shared_ptr< std::future< void > > > vcfWriterFutureFunctions;
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			auto vcfReaderPtr = iter.first;
			auto vcfPath = vcfReaderPtr->getFilePath();
			graphite::IFileWriter::SharedPtr fileWriter = outputFileMap[vcfPath];
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

    if (isIGVOutput)
    {
        // Write out GraphPaths fasta.
        std::string fastaFileName = firstFileName_withoutExtension + "." + "fa";
        updateFileMap(outputFileMap, fastaFileName, asciiFileType, outputDirectory, "fa");

        for (auto& hs: headerSequenceMap)
        {
            // Need +2 to account for the ">" and '\n' symbols.
            outputFileMap.at(fastaFileName)->write((">" + hs.first + '\n').c_str(), hs.first.length() + 2);
            outputFileMap.at(fastaFileName)->write((hs.second + '\n').c_str(), hs.second.length() + 1);
        }

        // Write out GraphPaths bed.
        std::string bedFileName = firstFileName_withoutExtension + "." + "bed";
        updateFileMap(outputFileMap, bedFileName, asciiFileType, outputDirectory, "bed");
        for (auto& hs : headerSequenceMap)
        {
            uint8_t nodeStringStart= hs.first.find("_") + 1;
            std::string nodeString = hs.first.substr(nodeStringStart);
            std::vector< std::string > nodeVector;
            if (nodeString.find(":") == std::string::npos)
                nodeVector.push_back(nodeString);
            else
                graphite::split(nodeString, ':', nodeVector);

            int32_t startPosition = 0;
            int32_t endPosition;
            std::string variantType;

            for (auto& node : nodeVector)
            {
                auto iter = nodeInfoMap.find(std::stoi(node));
                endPosition = startPosition + iter->second->getLength();
                if (iter->second->getVariantType() == 0)
                    variantType = "REF";
                else
                    variantType = "ALT";

                std::string bedLine = hs.first + '\t' + std::to_string(startPosition) + '\t' + std::to_string(endPosition) + '\t' + variantType + '\n';
                outputFileMap.at(bedFileName)->write(bedLine.c_str(), bedLine.length());

                startPosition = endPosition;
            }
        }
        
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

        std::string samFileName = firstFileName_withoutExtension + "." + "sam";
        updateFileMap(outputFileMap, samFileName, samFileType, outputDirectory, "sam");
        
        // Write out SAM header.
        graphite::SAMFileWriter::SharedPtr samFilePtr = std::dynamic_pointer_cast< graphite::SAMFileWriter >(outputFileMap.at(samFileName));
        samFilePtr->writeSamHeader(samHeader.c_str(), samHeader.length());

        // Adjust stream positions for the temporary SAM file and the final sam file.
        tempSamFilePtr->setInStreamToBeginning();

        // Write alignments from temp SAM file to final SAM file.
        while (!tempSamFilePtr->endOfFile())
        {
            std::string line;
            tempSamFilePtr->getInLine(line);
            samFilePtr->write(line.c_str(), line.length());
        }
    }
    
    // Close all open files.
	for (auto& iter : outputFileMap)
	{
		graphite::IFileWriter::SharedPtr fileWriter = iter.second;
		fileWriter->close();
	}

    // Remove temporary file containing alignments.
    if (isIGVOutput)
    {
        remove(tempSamFilePtr->getFilePath().c_str());
    }

	return 0;
}
