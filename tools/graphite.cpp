//#include <Python.h>
#include "core/util/Params.h"
#include "core/region/Region.h"
#include "core/reference/FastaReference.h"
#include "core/vcf/VCFReader.h"
#include "core/vcf/VCFWriter.h"
#include "core/bam/BamReader.h"
#include "core/graph/GraphProcessor.h"

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
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto includeDuplicates = params.getIncludeDuplicates();
	auto outputVisualizationFiles = params.outputVisualizationFiles();
	auto mappingQuality = params.getMappingQualityFilter();
	auto readSampleLimit = params.getReadSampleNumber();
	auto overwriteSampleName = params.getOverwrittenSampleName();

    // create reference reader
	auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath);

	// track samples from bams
	std::vector< graphite::Sample::SharedPtr > bamSamplePtrs;

	// create bam readers
	std::vector< graphite::BamReader::SharedPtr > bamReaderPtrs;
	for (auto bamPath : bamPaths)
	{
		auto bamReaderPtr = std::make_shared< graphite::BamReader >(bamPath);
		if (overwriteSampleName.length() == 0)
		{
			auto samplePtrs = bamReaderPtr->getSamplePtrs();
			bamSamplePtrs.insert(bamSamplePtrs.begin(), samplePtrs.begin(), samplePtrs.end());
		}
		else
		{
			auto samplePtr = std::make_shared< graphite::Sample >(overwriteSampleName, "1", bamPath);
			bamSamplePtrs.emplace_back(samplePtr);
			bamReaderPtr->overwriteSample(samplePtr);
		}
		bamReaderPtrs.emplace_back(bamReaderPtr);
	}

	// create VCF readers
	// create VCF writers
	std::vector< graphite::VCFReader::SharedPtr > vcfReaderPtrs;
	for (auto vcfPath : vcfPaths)
	{
		auto vcfWriterPtr = std::make_shared< graphite::VCFWriter >(vcfPath, bamSamplePtrs, outputDirectory);
		auto vcfReaderPtr = std::make_shared< graphite::VCFReader >(vcfPath, bamSamplePtrs, paramRegionPtr, vcfWriterPtr);
		vcfReaderPtrs.emplace_back(vcfReaderPtr);
	}

	// create graph processor
	// call process on processor
	auto graphProcessorPtr = std::make_shared< graphite::GraphProcessor >(fastaReferencePtr, bamReaderPtrs, vcfReaderPtrs, matchValue, misMatchValue, gapOpenValue, gapExtensionValue, outputVisualizationFiles, mappingQuality, readSampleLimit);
	graphProcessorPtr->processVariants();

	return 0;
}
