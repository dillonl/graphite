#include "core/alignments/BamAlignmentReader.h"
#include "core/alignments/BamAlignmentReaderManager.h"
#include "core/alignments/BamAlignmentReaderPreloadManager.h"
#include "core/variants/VCFFileReader.h"
#include "core/variants/IVariant.h"
#include "core/variants/VariantListVCFPreloaded.h"
#include "core/reference/FastaReference.h"
#include "core/utils/Parameters.h"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"
#include "core/utils/ThreadPool.hpp"

#include <thread>

int main(int argc, char** argv)
{
	gwiz::Parameters::Instance()->setParams(argc, argv);
	std::string fastaPath = gwiz::Parameters::Instance()->getFastaPath();
	std::string vcfPath = gwiz::Parameters::Instance()->getInVCFPath();
	std::string bamPath = gwiz::Parameters::Instance()->getBAMPath();
	std::string outputVCFPath = gwiz::Parameters::Instance()->getOutVCFPath();
	std::string region = gwiz::Parameters::Instance()->getRegion();
	if (gwiz::Parameters::Instance()->getThreadCount() > 0)
	{
		gwiz::ThreadPool::Instance()->setThreadCount(gwiz::Parameters::Instance()->getThreadCount());
	}

	gwiz::Region::SharedPtr regionPtr = std::make_shared< gwiz::Region >(region);
	auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);

	auto vcfFileReaderPtr = std::make_shared< gwiz::VariantListVCFPreloaded >(vcfPath, regionPtr);
	auto bamAlignmentReaderPreloadManager = std::make_shared< gwiz::BamAlignmentReaderPreloadManager >(bamPath, regionPtr);

	std::thread vcfLoadThread(&gwiz::VariantListVCFPreloaded::loadVariantsFromFile, vcfFileReaderPtr);
	std::thread loadBamsThread(&gwiz::BamAlignmentReaderPreloadManager::processBam, bamAlignmentReaderPreloadManager);
	vcfLoadThread.join();
	std::cout << "Finished loading vcf" << std::endl;
	loadBamsThread.join();

	std::cout << "Finished loading BAM and VCF" << std::endl;

	auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >();
	auto gsswGraphManager = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, vcfFileReaderPtr, bamAlignmentReaderPreloadManager, gsswAdjudicator);
	auto variantList = gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

	// std::cout << "starting to print vcf" << std::endl;

	if (outputVCFPath.size() > 0)
	{
		std::ofstream outVCF;
		outVCF.open(outputVCFPath, std::ios::out | std::ios::trunc);
		variantList->printToVCF(outVCF);
		outVCF.close();
	}
	else
	{
		variantList->printToVCF(std::cout);
	}
	return 0;
}
