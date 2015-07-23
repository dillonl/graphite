#include "core/alignment/BamAlignmentManager.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/util/Params.h"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"
#include "core/util/ThreadPool.hpp"

#include <thread>

int main(int argc, char** argv)
{
	gwiz::Params params;
	params.parseGSSW(argc, argv);
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

	// load bam alignments
	auto bamAlignmentManager = std::make_shared< gwiz::BamAlignmentManager >(bamPath, regionPtr);
	bamAlignmentManager->asyncLoadAlignments(); // begin the process of loading the alignments asynchronously

	// load variants from vcf
	auto variantManagerPtr = std::make_shared< gwiz::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously

	bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

	std::cout << "loaded vcfs and bams" << std::endl;

	// create an adjudicator for the graph
	auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >(swPercent);
	// the gsswGraphManager adjudicates on the variantManager's variants
	auto gsswGraphManager = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
	gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

	// get the complete variants list out of the variantListManager. The graphManager has adjudicated these variants.
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
