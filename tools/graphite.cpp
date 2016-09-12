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
#include <fstream>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "tabix.h"
#include "bgzf.h"

#include <zlib.h>

inline bool file_exists (const std::string& name)
{
	ifstream f(name.c_str());
	return f;
}
// void writeVariantListToFile(const std::string& path, graphite::VariantList::SharedPtr variantListPtr);
void validatePath(const std::string& path, const std::string& errorMessage, bool exitOnFailure)
{
	if (!file_exists(path))
	{
		std::cout << errorMessage << std::endl;
		if (exitOnFailure) {
			exit(EXIT_FAILURE);
		}
	}
}

void writeVariantListToFile(std::string path, graphite::VCFHeader::SharedPtr vcfHeaderPtr, graphite::VariantList::SharedPtr variantListPtr)
{
	bool fileExists = file_exists(path);
	std::ofstream outVCF;
	if (!fileExists)
	{
		outVCF.open(path, std::ios::out | std::ios::trunc);
	}
	else
	{
		outVCF.open(path, std::ios::out | std::ios::app);
	}

	variantListPtr->printToVCF(vcfHeaderPtr, !fileExists, outVCF);
	outVCF.close();
}

void writeVariantListToCompressedFile(std::string path, graphite::VCFHeader::SharedPtr vcfHeaderPtr, graphite::VariantList::SharedPtr variantListPtr)
{
	int fd = 0;


	if (file_exists(path))
	{
		fd = open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0666);
		variantListPtr->printToCompressedVCF(vcfHeaderPtr, false, fd);
	}
	else
	{
		fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
		variantListPtr->printToCompressedVCF(vcfHeaderPtr, true, fd);
	}
}

int main(int argc, char** argv)
{
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
	if (outputDirectory.size() == 0)
	{
		validatePath(outputDirectory, "Invalid output directory, please create that directory or change to an existing directory and rerun Graphite", true);
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
		std::string path = vcfPath.substr(vcfPath.find_last_of("/") + 1);
		std::string filePath = vcfoutPaths[vcfPath] = outputDirectory + "/" + path;
		uint32_t counter = 1;
		while (file_exists(filePath + ".gz"))
		{
			std::string extension = vcfPath.substr(vcfPath.find_last_of(".") + 1);
			std::string fileNameWithoutExtension = path.substr(0, path.find_last_of("."));
			filePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + extension;
			++counter;
		}
		filePath += ".tmp";
		vcfoutPaths[vcfPath] = filePath;
	}

	std::unordered_set< std::string > outputPaths;
	bool firstTime = true;

	for (uint32_t regionCount = 0; regionCount < regionPtrs.size(); ++regionCount)
	{
		auto regionPtr = regionPtrs[regionCount];

		auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

		// load variants from vcf
		auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr, fastaReferencePtr, graphSize);
		variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously

		variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
		// std::cout << "loaded vcf region: " << regionPtr->getRegionString() << std::endl;
		// load bam alignments
		auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(samplePtrs, regionPtr, excludeDuplicates);
		bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
		bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory

		// std::cout << "loaded alignments region: " << regionPtr->getRegionString() << std::endl;

		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
		bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

		milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

		// create an adjudicator for the graph
		auto gsswAdjudicator = std::make_shared< graphite::adjudicator::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);

		// the gsswGraphManager adjudicates on the variantManager's variants
		auto gsswGraphManager = std::make_shared< graphite::adjudicator::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), graphSize, 1000, 100);

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

		graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);

		std::vector< std::shared_ptr< std::thread > > fileWriters;
		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFReadersAndVariantListsMap();
		std::deque< std::shared_ptr< std::future< void > > > vcfWriterFutureFunctions;
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			auto vcfReaderPtr = iter.first;
			auto vcfPath = vcfReaderPtr->getFilePath();
			std::string currentVCFOutPath = vcfoutPaths[vcfPath];
			auto variantListPtr = iter.second;
			auto vcfHeaderPtr = vcfReaderPtr->getVCFHeader();
			for (auto samplePtr : bamAlignmentManager->getSamplePtrs())
			{
				vcfHeaderPtr->registerSample(samplePtr);
			}
			if (firstTime)
			{
				outputPaths.emplace(currentVCFOutPath);
			}
			// auto funct = std::bind(&writeVariantListToCompressedFile, currentVCFOutPath, vcfHeaderPtr, variantListPtr);
			auto funct = std::bind(&writeVariantListToFile, currentVCFOutPath, vcfHeaderPtr, variantListPtr);
			auto functFuture = graphite::ThreadPool::Instance()->enqueue(funct);
			vcfWriterFutureFunctions.push_back(functFuture);
		}

		while (!vcfWriterFutureFunctions.empty())
		{
			vcfWriterFutureFunctions.front()->wait();
			vcfWriterFutureFunctions.pop_front();
		}

		graphite::MappingManager::Instance()->clearRegisteredMappings();
		firstTime = false;
	}

	for (auto& outputPath : outputPaths)
	{
		std::string tmpOutputPath = outputPath.substr(0, outputPath.find_last_of(".")); // remove the .tmp
		std::string extension = tmpOutputPath.substr(tmpOutputPath.find_last_of(".") + 1); // get the extension
		if (extension.compare("gz") == 0) // if it's a gz extension then remove it
		{
			tmpOutputPath = outputPath.substr(0, tmpOutputPath.find_last_of("."));
			extension = tmpOutputPath.substr(tmpOutputPath.find_last_of(".") + 1);
		}
		std::string destOutputPath = tmpOutputPath += ".gz";
		int fp;
		BGZF* f_dst;
		void* buffer;
		long start = 0;
		long end = -1;
		int c;
		int is_forced = 0;
		int WINDOW_SIZE = 64 * 1024;
		fp = open(outputPath.c_str(), O_RDONLY);
		f_dst = bgzf_open(destOutputPath.c_str(), "w");
		buffer = malloc(WINDOW_SIZE);
		while ((c = read(fp, buffer, WINDOW_SIZE)) > 0)
		{
			if (bgzf_write(f_dst, buffer, c) < 0)
			{
				std::cout << "an error occurred while compressing" << std::endl;
				exit(0);
			}
		}
		free(buffer);
		close(fp);

		bgzf_close(f_dst);

		ti_conf_t* conf = &ti_conf_vcf;
		ti_index_build(destOutputPath.c_str(), conf);
		std::string indexPath = destOutputPath + ".tbi";
		std::string tmpIndexPath = destOutputPath.substr(0, outputPath.find_last_of(".")) + ".tbi"; // rename tbi
		rename(indexPath.c_str(), tmpIndexPath.c_str());
		remove(outputPath.c_str());
	}
	return 0;
}
