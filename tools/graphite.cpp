#include "core/alignment/AlignmentManager.hpp"
#include "core/alignment/HTSLibAlignmentReader.h"
#include "core/alignment/HTSLibAlignmentManager.h"
// #include "core/alignment/BamAlignmentManager.h"
// #include "core/alignment/BamAlignmentReader.h"
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
#include "core/file/BGZFFileWriter.h"
#include "core/file/ASCIIFileWriter.h"

#include <thread>
#include <unordered_set>
#include <fstream>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

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
	// graphite::AlignmentManager< HTSLibAlignmentReader > tmp;
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
	graphite::FileType fileType = graphite::FileType::ASCII;

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
		// auto tmpSamplePtrs = graphite::BamAlignmentReader::GetBamReaderSamples(bamPath);
		auto tmpSamplePtrs = graphite::HTSLibAlignmentReader::GetBamReaderSamples(bamPath);
		samplePtrs.insert(samplePtrs.end(), tmpSamplePtrs.begin(), tmpSamplePtrs.end());
	}
	for (auto samplePtr : samplePtrs)
	{
		graphite::SampleManager::Instance()->addSamplePtr(samplePtr);
	}

	auto alignmentReaderManagerPtr = std::make_shared< graphite::AlignmentReaderManager< graphite::HTSLibAlignmentReader > >(bamPaths, threadCount);

	// make sure paths exist
	validatePath(fastaPath, "Invalid Fasta path, please provide the correct path to the Fasta file and rerun Graphite", true);

	std::unordered_map< std::string, graphite::IFileWriter::SharedPtr > vcfoutPaths;
	for (auto vcfPath : vcfPaths)
	{
		std::string path = vcfPath.substr(vcfPath.find_last_of("/") + 1);
		std::string filePath = outputDirectory + "/" + path;
		uint32_t counter = 1;
		while (file_exists(filePath))
		{
			std::string extension = vcfPath.substr(vcfPath.find_last_of(".") + 1);
			std::string fileNameWithoutExtension = path.substr(0, path.find_last_of("."));
			filePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + extension;
			++counter;
		}
		// filePath += ".tmp";
		graphite::IFileWriter::SharedPtr fileWriterPtr;
		if (fileType == graphite::FileType::BGZF)
		{
			fileWriterPtr = std::make_shared< graphite::BGZFFileWriter >(filePath);
		}
		else
		{
			fileWriterPtr = std::make_shared< graphite::ASCIIFileWriter >(filePath);
		}
		fileWriterPtr->open();
		vcfoutPaths[vcfPath] = fileWriterPtr;
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
		auto alignmentManager = std::make_shared< graphite::HTSLibAlignmentManager >(samplePtrs, regionPtr, alignmentReaderManagerPtr, excludeDuplicates);
		alignmentManager->loadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
		// auto bamAlignmentManager = std::make_shared< graphite::BamAlignmentManager >(samplePtrs, regionPtr, excludeDuplicates);
		// bamAlignmentManager->asyncLoadAlignments(variantManagerPtr, graphSize); // begin the process of loading the alignments asynchronously
		// bamAlignmentManager->waitForAlignmentsToLoad(); // wait for alignments to load into memory
		// std::cout << "loaded alignments region: " << regionPtr->getRegionString() << std::endl;

		variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
		// bamAlignmentManager->releaseResources(); // release the bam file into memory, we no longer need the file resources

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

		milliseconds_since_epoch = std::chrono::system_clock::now().time_since_epoch() / std::chrono::milliseconds(1);

		// create an adjudicator for the graph
		auto gsswAdjudicator = std::make_shared< graphite::adjudicator::GSSWAdjudicator >(swPercent, matchValue, misMatchValue, gapOpenValue, gapExtensionValue);

		// the gsswGraphManager adjudicates on the variantManager's variants
		// auto gsswGraphManager = std::make_shared< graphite::adjudicator::GraphManager >(fastaReferencePtr, variantManagerPtr, bamAlignmentManager, gsswAdjudicator);
		auto gsswGraphManager = std::make_shared< graphite::adjudicator::GraphManager >(fastaReferencePtr, variantManagerPtr, alignmentManager, gsswAdjudicator);
		gsswGraphManager->buildGraphs(fastaReferencePtr->getRegion(), graphSize, 1000, 100);

		graphite::MappingManager::Instance()->evaluateAlignmentMappings(gsswAdjudicator);
		graphite::MappingManager::Instance()->clearRegisteredMappings();

		std::vector< std::shared_ptr< std::thread > > fileWriters;
		auto vcfPathsAndVariantListPtrsMap = variantManagerPtr->getVCFReadersAndVariantListsMap();
		std::deque< std::shared_ptr< std::future< void > > > vcfWriterFutureFunctions;
		for (auto& iter : vcfPathsAndVariantListPtrsMap)
		{
			auto vcfReaderPtr = iter.first;
			auto vcfPath = vcfReaderPtr->getFilePath();
			graphite::IFileWriter::SharedPtr fileWriter = vcfoutPaths[vcfPath];
			std::string currentVCFOutPath = fileWriter->getFilePath();
			auto variantListPtr = iter.second;
			auto vcfHeaderPtr = vcfReaderPtr->getVCFHeader();

            // for (auto samplePtr : bamAlignmentManager->getSamplePtrs())
			// for (auto samplePtr : alignmentManager->getSamplePtrs())
			for (auto samplePtr : samplePtrs)
			{
				vcfHeaderPtr->registerSample(samplePtr);
			}
			if (firstTime)
			{
				outputPaths.emplace(currentVCFOutPath);
			}
			auto funct = std::bind(&graphite::VariantList::writeVariantList, variantListPtr, fileWriter, vcfHeaderPtr, firstTime);
			auto functFuture = graphite::ThreadPool::Instance()->enqueue(funct);
			vcfWriterFutureFunctions.push_back(functFuture);
		}

		while (!vcfWriterFutureFunctions.empty())
		{
			vcfWriterFutureFunctions.front()->wait();
			vcfWriterFutureFunctions.pop_front();
		}

		graphite::SequenceManager::Instance()->clearSequences();
		firstTime = false;
	}

	for (auto& iter : vcfoutPaths)
	{
		graphite::IFileWriter::SharedPtr fileWriter = iter.second;
		fileWriter->close();
	}

	graphite::adjudicator::GSSWAdjudicator* adj_p;
	std::cout << "adj counts: " << (uint32_t)adj_p->s_adj_count << " [total]" << std::endl;

	return 0;
}
