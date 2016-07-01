#include "core/alignment/BamAlignmentManager.h"
#include "core/alignment/SamtoolsAlignmentReader.h"
#include "core/alignment/BamAlignmentReader.h"
#include "core/variant/VCFManager.h"
#include "core/variant/VCFFileReader.h"
#include "core/reference/FastaReference.h"
#include "core/mapping/MappingManager.h"
#include "core/variant/VCFHeader.h"
#include "core/util/Params.h"
#include "core/util/ThreadPool.hpp"
#include "adjudicator/graph/GraphManager.h"
#include "adjudicator/graph/GSSWAdjudicator.h"
#include "core/variant/VCFHeader.h"
#include "core/alignment/SampleManager.hpp"
#include "core/region/Region.h"

#include <thread>
#include <unordered_set>

#include <boost/filesystem.hpp>

// void writeVariantListToFile(const std::string& path, graphite::VariantList::SharedPtr variantListPtr);
void writeVariantListToFile(std::string path, bool writeHeader, graphite::VCFHeader::SharedPtr vcfHeaderPtr, graphite::VariantList::SharedPtr variantListPtr)
{
	std::ofstream outVCF;
	if (writeHeader)
	{
		outVCF.open(path, std::ios::out | std::ios::trunc);
	}
	else
	{
		outVCF.open(path, std::ios::out | std::ios::app);
	}

	variantListPtr->printToVCF(vcfHeaderPtr, writeHeader, outVCF);
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
	std::vector< std::string > inclusionList = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"};
	unsigned long milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);
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
	auto paramRegionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto excludeDuplicates = params.getExcludeDuplicates();
	auto graphSize = params.getGraphSize();
	graphite::ThreadPool::Instance()->setThreadCount(threadCount);

	std::vector< graphite::Region::SharedPtr > regionPtrs;
	if (paramRegionPtr == nullptr)
	{
		std::unordered_set< std::string > regionStringSet;
		for (auto vcfPath : vcfPaths)
		{
			auto tmpRegionPtrs = graphite::VCFFileReader::GetAllRegionsInVCF(vcfPath);
			for (auto& tmpRegionPtr : tmpRegionPtrs)
			{
				if (regionStringSet.find(tmpRegionPtr->getRegionString()) == regionStringSet.end())
				{
					regionStringSet.emplace(tmpRegionPtr->getRegionString());
					regionPtrs.emplace_back(tmpRegionPtr);
				}
			}
		}
	}
	else
	{
		regionPtrs.emplace_back(paramRegionPtr);
	}

	// skip the bampath checking for now
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

	std::vector< graphite::Sample::SharedPtr > samplePtrs;
	for (auto bamPath : bamPaths)
	{
		// auto tmpSamplePtrs = graphite::SamtoolsAlignmentReader::GetBamReaderSamples(bamPath);
		auto tmpSamplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
		samplePtrs.insert(samplePtrs.end(), tmpSamplePtrs.begin(), tmpSamplePtrs.end());
	}
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}

	// make sure paths exist
	validatePath(fastaPath, "Invalid Fasta path, please provide the correct path to the Fasta file and rerun Graphite", true);

	std::unordered_map< std::string, std::string > vcfoutPaths;
	for (auto vcfPath : vcfPaths)
	{

		boost::filesystem::path vcfFSPath(vcfPath);
		boost::filesystem::path outputDirectoryPath(outputDirectory);
		boost::filesystem::path extension(".vcf");
		boost::filesystem::path outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfFSPath.stem().string() + extension.string());
		uint32_t counter = 1;
		while (boost::filesystem::exists(outputVCFPath.string()))
		{
			boost::filesystem::path countPath("_" + std::to_string(counter++));
			outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfFSPath.stem().string() + countPath.string() + extension.string());
		}
		vcfoutPaths[vcfPath] = outputVCFPath.string();
	}

	for (uint32_t regionCount = 0; regionCount < regionPtrs.size(); ++regionCount)
	{
		auto regionPtr = regionPtrs[regionCount];
		if (std::find(inclusionList.begin(), inclusionList.end(), regionPtr->getReferenceID()) == inclusionList.end())
		{
			std::cout << "skipping: " << regionPtr->getRegionString() << std::endl;
			continue;
		}
		std::cout << "processing: " << regionPtr->getRegionString() << std::endl;

		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

		// load variants from vcf
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, graphSize);
		variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously

		variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
		// load bam alignments
		auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(samplePtrs, regionPtr, excludeDuplicates);
		bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
		bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory
		std::cout << "finished reading alignments" << std::endl;

		for (auto samplePtr : bamAlignmentManager->getSamplePtrs())
		{
			vcfHeaderPtr->registerSample(samplePtr);
		}

		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
		bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

		milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

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
			variantListPtr->printToVCF(vcfHeaderPtr, (regionCount == 0), std::cout);
		}
		else
		{
			std::vector< std::shared_ptr< std::thread > > fileWriters;
			auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFPathsAndVariantListsMap();
			for (auto& iter : vcfPathsAndVariantListPtrsMap)
			{
				boost::filesystem::path vcfPath(iter.first);
				std::string currentVCFOutPath = vcfoutPaths[iter.first];
				auto variantListPtr = iter.second;
				/*
				boost::filesystem::path outputDirectoryPath(outputDirectory);
				boost::filesystem::path extension(".vcf");
				boost::filesystem::path outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfPath.stem().string() + extension.string());
				uint32_t counter = 1;
				while (boost::filesystem::exists(outputVCFPath.string()))
				{
					boost::filesystem::path countPath("_" + std::to_string(counter++));
					outputVCFPath = outputDirectoryPath / boost::filesystem::path(vcfPath.stem().string() + countPath.string() + extension.string());
				}
				*/
				// auto funct = std::bind(&writeVariantListToFile, outputVCFPath.string(), (regionCount == 0), vcfHeaderPtr, variantListPtr);
				auto funct = std::bind(&writeVariantListToFile, currentVCFOutPath, (regionCount == 0), vcfHeaderPtr, variantListPtr);
				auto functFuture = graphite::ThreadPool::Instance()->enqueue(funct);
			}

			graphite::ThreadPool::Instance()->joinAll();
		}
		graphite::MappingManager::Instance()->clearRegisteredMappings();
	}
	return 0;
}
