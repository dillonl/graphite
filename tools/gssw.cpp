#include "core/alignment/BamAlignmentReader.h"
#include "core/alignment/BamAlignmentReaderManager.h"
#include "core/alignment/BamAlignmentReaderPreloadManager.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/util/Parameters.h"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"
#include "core/util/ThreadPool.hpp"

#include <thread>

int main(int argc, char** argv)
{
	gwiz::Parameters::Instance()->setParams(argc, argv);
	auto fastaPath = gwiz::Parameters::Instance()->getFastaPath();
	auto vcfPaths = gwiz::Parameters::Instance()->getInVCFPaths();
	auto bamPath = gwiz::Parameters::Instance()->getBAMPath();
	auto outputVCFPath = gwiz::Parameters::Instance()->getOutVCFPath();
	auto region = gwiz::Parameters::Instance()->getRegion();
	if (gwiz::Parameters::Instance()->getThreadCount() > 0)
	{
		gwiz::ThreadPool::Instance()->setThreadCount(gwiz::Parameters::Instance()->getThreadCount());
	}

	auto regionPtr = std::make_shared< gwiz::Region >(region);
	auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);

	auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);

	std::thread loadBamsThread(&gwiz::BamAlignmentReaderPreloadManager::processBam, bamAlignmentReaderPreloadManager);
	auto variantManagerPtr = std::make_shared< gwiz::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs();
	loadBamsThread.join();
	variantManagerPtr->waitForVCFsToLoadAndProcess();

	auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >(gwiz::Parameters::Instance()->getSWPercent());
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
