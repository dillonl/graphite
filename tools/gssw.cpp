#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/mapping/MappingManager.h"
#include "core/util/Params.h"
#include "core/util/ThreadPool.hpp"
#include "gssw/graph/GraphManager.h"
#include "gssw/graph/GSSWAdjudicator.h"

#include <thread>

#include <boost/filesystem.hpp>

// void writeVariantListToFile(const std::string& path, gwiz::VariantList::SharedPtr variantListPtr);
void writeVariantListToFile(std::string path, gwiz::VariantList::SharedPtr variantListPtr)
{
	std::ofstream outVCF;
	outVCF.open(path, std::ios::out | std::ios::trunc);
	variantListPtr->printToVCF(outVCF);
	outVCF.close();
}

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
	auto outputDirectory = params.getOutputDirectory();
	auto regionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
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

	// create an adjudicator for the graph
	auto gsswAdjudicator = std::make_shared< gwiz::gssw::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);
	// the gsswGraphManager adjudicates on the variantManager's variants
	auto gsswGraphManager = std::make_shared< gwiz::gssw::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
	gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), 3000, 1000, 100);

	gwiz::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);

	// get the complete variants list out of the variantListManager. The graphManager has adjudicated these variants.

	gwiz::ThreadPool::Instance()->joinAll();

	if (outputDirectory.size() == 0)
	{
		auto variantListPtr = variantManagerPtr->getCompleteVariantList();
		variantListPtr->printToVCF(std::cout);
	}
	else
	{
		std::vector< std::shared_ptr< std::thread > > fileWriters;
		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFPathsAndVariantListsMap();
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			boost::filesystem::path vcfPath(iter.first);
			auto variantListPtr = iter.second;
			boost::filesystem::path outputDirectoryPath(outputDirectory);
			boost::filesystem::path extension(".vcf");
			boost::filesystem::path outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfPath.stem().string() + extension.string());
			uint32_t counter = 1;
			while (boost::filesystem::exists(outputVCFPath.string()))
			{
				boost::filesystem::path countPath("_" + std::to_string(counter++));
				outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfPath.stem().string() + countPath.string() + extension.string());
			}
			auto fileWriterThread = std::make_shared< std::thread >(writeVariantListToFile, outputVCFPath.string(), variantListPtr);
			fileWriters.emplace_back(fileWriterThread);
		}
		for (auto& fileWriterThread : fileWriters)
		{
			fileWriterThread->join();
		}
		fileWriters.clear(); // clear file writer threads
	}
	return 0;
}
