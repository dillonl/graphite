#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/mapping/MappingManager.h"
#include "core/variant/VCFHeader.h"
#include "core/util/Params.h"
#include "core/util/ThreadPool.hpp"
#include "adjudicator/graph/GraphManager.h"
#include "adjudicator/graph/GSSWAdjudicator.h"
#include "core/variant/VCFHeader.h"
#include "core/alignment/SampleManager.hpp"

#include <thread>

#include <boost/filesystem.hpp>

// void writeVariantListToFile(const std::string& path, graphite::VariantList::SharedPtr variantListPtr);
void writeVariantListToFile(std::string path, graphite::VCFHeader::SharedPtr vcfHeaderPtr, graphite::VariantList::SharedPtr variantListPtr)
{
	std::ofstream outVCF;
	outVCF.open(path, std::ios::out | std::ios::trunc);

	variantListPtr->printToVCF(vcfHeaderPtr, outVCF);
	outVCF.close();
}

void validatePath(const std::string& path, const std::string& errorMessage, bool exitOnFailure)
{
	if (!boost::filesystem::exists(path))
	{
		std::cout << errorMessage << std::endl;
		if (exitOnFailure) {
			exit(EXIT_FAILURE);
		}
	}
}

int main(int argc, char** argv)
{
	unsigned long milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << "starting: " << milliseconds_since_epoch << std::endl;
	graphite::Params params;
	params.parseGSSW(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		exit(0);
	}
	auto fastaPath = params.getFastaPath();
	auto vcfPaths = params.getInVCFPaths();
	auto bamPaths = params.getBAMPaths();
	auto outputDirectory = params.getOutputDirectory();
	auto regionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto graphSize = params.getGraphSize();
	graphite::ThreadPool::Instance()->setThreadCount(threadCount);

	// make sure paths exist
	validatePath(fastaPath, "Invalid Fasta path, please provide the correct path to the Fasta file and rerun Graphite", true);
	for (auto bamPath : bamPaths)
	{
		validatePath(bamPath, "Invalid BAM path: " + bamPath + ", please provide the correct path to the BAM and rerun Graphite", true);
	}
	for (auto vcfPath : vcfPaths)
	{
		validatePath(vcfPath, "Invalid VCF path: " + vcfPath + ", please provide the correct path to the VCF and rerun Graphite", true);
	}
	if (outputDirectory.size() > 0)
	{
		validatePath(outputDirectory, "Invalid output directory, please create that directory or change to an existing directory and rerun Graphite", true);
	}

	graphite::VCFHeader::SharedPtr vcfHeaderPtr = std::make_shared< graphite::VCFHeader >();

	for (auto bamPath : bamPaths)
	{
		vcfHeaderPtr->addHeaderLine("##samplefile=" + bamPath);
	}

	auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

	std::vector< graphite::Sample::SharedPtr > samplePtrs;
	for (auto bamPath : bamPaths)
	{
		auto tmpSamplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
		samplePtrs.insert(samplePtrs.end(), tmpSamplePtrs.begin(), tmpSamplePtrs.end());
	}
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}

	// load variants from vcf
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously

	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	// load bam alignments
	auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(samplePtrs, regionPtr);
	bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
	bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory

	for (auto samplePtr : bamAlignmentManager->getSamplePtrs())
	{
		vcfHeaderPtr->registerSample(samplePtr);
	}

	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

	milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
	std::cout << "ending: " << milliseconds_since_epoch << std::endl;

	// create an adjudicator for the graph
	auto gsswAdjudicator = std::make_shared< graphite::adjudicator::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);
	// the gsswGraphManager adjudicates on the variantManager's variants
	auto gsswGraphManager = std::make_shared< graphite::adjudicator::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
	gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), graphSize, 1000, 100);

	if (outputDirectory.size() == 0)
	{
		variantManagerPtr->getCompleteVariantList()->processOverlappingAlleles();
	}
	else
	{
		for (auto& iter : variantManagerPtr->getVCFPathsAndVariantListsMap())
		{
			graphite::ThreadPool::Instance()->enqueue(std::bind(&graphite::IVariantList::processOverlappingAlleles, iter.second));
		}
	}
	graphite::ThreadPool::Instance()->joinAll();

	graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);

	// get the complete variants list out of the variantListManager. The graphManager has adjudicated these variants.

	graphite::ThreadPool::Instance()->joinAll();

	if (outputDirectory.size() == 0)
	{
		auto variantListPtr = variantManagerPtr->getCompleteVariantList();
		variantListPtr->printToVCF(vcfHeaderPtr, std::cout);
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
			auto funct = std::bind(&writeVariantListToFile, outputVCFPath.string(), vcfHeaderPtr, variantListPtr);
			auto functFuture = graphite::ThreadPool::Instance()->enqueue(funct);
		}

		graphite::ThreadPool::Instance()->joinAll();
	}
	return 0;
}
