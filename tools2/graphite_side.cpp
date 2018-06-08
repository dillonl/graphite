//#include <Python.h>
#include "core2/util/Params.h"
#include "core2/region/Region.h"
#include "core2/reference/FastaReference.h"
#include "core2/vcf/VCFReader.h"
#include "core2/vcf/VCFWriter.h"
#include "core2/bam/BamReader.h"
#include "core2/graph/GraphProcessor.h"

#include <string>
#include <iostream>

int main(int argc, char** argv)
{
	// get all param properties
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
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto includeDuplicates = params.getIncludeDuplicates();
	auto outputVisualizationFiles = params.outputVisualizationFiles();

    // create reference reader
	auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath);

	// track samples from bams
	std::vector< graphite::Sample::SharedPtr > bamSamplePtrs;

	// create bam readers
	std::vector< graphite::BamReader::SharedPtr > bamReaderPtrs;
	for (auto bamPath : bamPaths)
	{
		auto bamReaderPtr = std::make_shared< graphite::BamReader >(bamPath);
        auto samplePtrs = bamReaderPtr->getSamplePtrs();
		bamSamplePtrs.insert(bamSamplePtrs.begin(), samplePtrs.begin(), samplePtrs.end());
		bamReaderPtrs.emplace_back(bamReaderPtr);
	}

	// create VCF readers
	// create VCF writers
	std::vector< graphite::VCFReader::SharedPtr > vcfReaderPtrs;
	for (auto vcfPath : vcfPaths)
	{
		auto vcfReaderPtr = std::make_shared< graphite::VCFReader >(vcfPath, bamSamplePtrs, paramRegionPtr);
		auto vcfWriterPtr = std::make_shared< graphite::VCFWriter >(vcfPath, outputDirectory);
		vcfReaderPtr->registerVCFWriter(vcfWriterPtr);
		vcfReaderPtrs.emplace_back(vcfReaderPtr);
	}

	// create graph processor
	// call process on processor
	auto graphProcessorPtr = std::make_shared< graphite::GraphProcessor >(fastaReferencePtr, bamReaderPtrs, vcfReaderPtrs);
	graphProcessorPtr->processVariants();

	return 0;
}
