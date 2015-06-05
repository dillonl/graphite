#include "core/alignment/BamAlignmentReader.h"
#include "core/alignment/BamAlignmentReaderManager.h"
#include "core/alignment/BamAlignmentReaderPreloadManager.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/util/Params.h"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"
#include "core/util/ThreadPool.hpp"

#include <thread>

int main(int argc, char** argv)
{
	gwiz::Params params(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		exit(0);
	}
	auto fastaPath = params.getFastaPath();
	auto vcfPaths = params.getInVCFPaths();
	auto bamPath = params.getBAMPath();
	auto outputVCFPath = params.getOutVCFPath();
	auto regionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	gwiz::ThreadPool::Instance()->setThreadCount(threadCount);

	auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);
	auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);

	std::thread loadBamsThread(&gwiz::BamAlignmentReaderPreloadManager::processBam, bamAlignmentReaderPreloadManager);
	auto variantManagerPtr = std::make_shared< gwiz::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs();
	loadBamsThread.join();
	variantManagerPtr->waitForVCFsToLoadAndProcess();
	variantManagerPtr->releaseVCFResources(); // releases the vcf file memory

	std::cout << "loaded" << std::endl;

	auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >(swPercent);
	auto gsswGraphManager = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentReaderPreloadManager, gsswAdjudicator);
	gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

	auto variantListPtr = variantManagerPtr->getCompleteVariantList();
	if (outputVCFPath.size() > 0)
	{
		std::ofstream outVCF;
		outVCF.open(outputVCFPath, std::ios::out | std::ios::trunc);
		variantListPtr->printToVCF(outVCF);
		outVCF.close();
	}
	else
	{
		variantListPtr->printToVCF(std::cout);
	}
	return 0;
}
