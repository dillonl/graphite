#include "core/alignment/AlignmentManager.hpp"
#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/variant/VCFFileReader.h"
#include "core/reference/FastaReference.h"
#include "core/mapping/MappingManager.h"
#include "core/variant/VCFHeader.h"
#include "core/variant/VCFFileWriterManager.hpp"
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
#include "core/util/VisualizationToolKit.h"

#include <thread>
#include <unordered_set>
#include <fstream>
#include <stdio.h>

#include <zlib.h>

int main(int argc, char** argv)
{
	// graphite::AlignmentManager< HTSLibAlignmentReader > tmp;
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
	auto includeDuplicates = params.getIncludeDuplicates();
	auto graphSize = params.getGraphSize();
	auto outputVisualizationFiles = params.outputVisualizationFiles();
	graphite::FileType fileType = graphite::FileType::ASCII;

	// graphite::ThreadPool::Instance()->setThreadCount(threadCount);
	graphite::ThreadPool::Instance()->setThreadCount(1);

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
	graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(bamPaths);
	graphite::VCFFileWriterManager::Instance()->addVCFFileWritersForVCFs(vcfPaths, outputDirectory);

	graphite::VisualizationToolKit::SharedPtr vtkPtr = nullptr;
	if (outputVisualizationFiles)
	{
		vtkPtr = std::make_shared< graphite::VisualizationToolKit >(outputDirectory, bamPaths, matchValue);
	}

	for (uint32_t regionCount = 0; regionCount < regionPtrs.size(); ++regionCount)
	{
		// std::cout << "paused" << std::endl;
		// char tmp[256];
		// std::cin.getline(tmp, 256);

		// auto alignmentReaderManagerPtr = std::make_shared< graphite::AlignmentReaderManager< graphite::BamAlignmentReader > >(bamPaths, threadCount); // this used to go above this loop but it caused issues with loading bam regions from out-of-order VCFs
		auto regionPtr = regionPtrs[regionCount];
		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

		// load variants from vcf
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, readLength);
		variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
		variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources

		std::deque< std::shared_ptr< std::future< void > > > variantManagerFutureFunctions;
		for (auto& iter : variantManagerPtr->getVCFReadersAndVariantListsMap())
		{
			iter.second->processOverlappingAlleles();
		}

		// create an adjudicator for the graph
		auto gsswAdjudicator = std::make_shared< graphite::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);
		auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, bamPaths, sampleManagerPtr, false, includeDuplicates, gsswAdjudicator);
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength, vtkPtr);

		// graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);
		// graphite::MappingManager::Instance()->clearRegisteredMappings();

		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFReadersAndVariantListsMap();
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			graphite::VCFFileWriter::SharedPtr vcfFileWriterPtr = graphite::VCFFileWriterManager::Instance()->getVCFFileWriter(iter.first->getFilePath());
			auto vcfHeaderPtr = iter.first->getVCFHeader();
			vcfHeaderPtr->registerActiveSample(sampleManagerPtr);
			vcfFileWriterPtr->writeVariantList(iter.second, vcfHeaderPtr); // iter.second is the VariantList
		}
	}

	if (vtkPtr != nullptr)
	{
		vtkPtr->closeResources();
	}
	graphite::VCFFileWriterManager::Instance()->closeAllVCFFileWriters();

	return 0;
}
