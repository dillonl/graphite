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

#include <thread>
#include <unordered_set>
#include <fstream>
#include <stdio.h>

#include <zlib.h>

int main(int argc, char** argv)
{
	// graphite::AlignmentManager< HTSLibAlignmentReader > tmp;
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

	for (auto& iter : vcfoutPaths)
	{
		graphite::IFileWriter::SharedPtr fileWriter = iter.second;
		fileWriter->close();
	}

	// graphite::GSSWAdjudicator* adj_p;
	// std::cout << "adj counts: " << (uint32_t)adj_p->s_adj_count << " [total]" << std::endl;

	return 0;
}
