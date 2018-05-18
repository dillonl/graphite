#include "core/alignment/AlignmentManager.hpp"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/variant/VCFFileReader.h"
#include "core/reference/FastaReference.h"
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
	graphite::Params params;
	params.parseGSSW(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		return 0;
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

	// graphite::ThreadPool::Instance()->setThreadCount(1);

	std::vector< graphite::Region::SharedPtr > regionPtrs;
	if (paramRegionPtr == nullptr)
	{
		regionPtrs = graphite::VCFFileReader::GetAllRegionsInVCF(vcfPaths);
	}
	else
	{
		regionPtrs.emplace_back(paramRegionPtr);
	}

	uint32_t readLength = 0;
	for (auto bamPath : bamPaths)
	{
		uint32_t tmp = graphite::BamAlignmentReader::GetReadLength(bamPath);
		readLength = (tmp > readLength) ? tmp : readLength;
	}
	graphite::SampleManager::SharedPtr sampleManagerPtr = std::make_shared< graphite::SampleManager >(bamPaths);
	graphite::VCFFileWriterManager::Instance()->addVCFFileWritersForVCFs(vcfPaths, outputDirectory);

	graphite::VisualizationToolKit::SharedPtr vtkPtr = nullptr;
	if (outputVisualizationFiles)
	{
		vtkPtr = std::make_shared< graphite::VisualizationToolKit >(outputDirectory, bamPaths, matchValue);
	}

	for (uint32_t regionCount = 0; regionCount < regionPtrs.size(); ++regionCount)
	{
		auto regionPtr = regionPtrs[regionCount];
		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

		// load variants from vcf
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, readLength);
		variantManagerPtr->loadVariants(); // load the variants from the vcfs

		std::deque< std::shared_ptr< std::future< void > > > variantManagerFutureFunctions;
		for (auto& iter : variantManagerPtr->getVCFReadersAndVariantListsMap())
		{
			iter.second->processOverlappingAlleles();
		}

		// create an adjudicator for the graph
		auto gsswAdjudicator = std::make_shared< graphite::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);
		auto gsswGraphManager = std::make_shared< graphite::GraphManager >(fastaReferencePtr, variantManagerPtr, bamPaths, sampleManagerPtr, false, includeDuplicates, gsswAdjudicator);
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), readLength, vtkPtr);

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
