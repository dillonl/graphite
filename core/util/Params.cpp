#include "Params.h"
#include "Utility.h"
#include "config/GraphiteConfig.hpp"

#include <string.h>
#include <thread>
#include <iostream>

namespace graphite
{

	Params::Params() :
		m_options("Graphite", "VCF Allele Adjudicator")
	{

	}

	Params::~Params()
	{
	}

	void Params::parseGSSW(int argc, char** argv)
	{
		this->m_options.add_options()
			("h,help","Print help message")
			("d,include_duplicates", "Include Duplicate Reads")
			("v,vcf", "Path to input VCF file[s], separate multiple files by space", cxxopts::value< std::vector< std::string > >())
			("b,bam", "Path to input BAM file[s], separate multiple files by space", cxxopts::value< std::vector< std::string > >())
			("r,region", "Region information", cxxopts::value< std::string >())
			("o,output_directory", "Path to output directory", cxxopts::value< std::string >())
			("s,sample_name", "Ignore the BAM file's samples and use this passed in value as the sample name", cxxopts::value< std::string >()->default_value(""))
			("f,fasta", "Path to input FASTA file", cxxopts::value< std::string >())
			("p,percent_match", "Smith-Waterman Percent [optional - default is 90]", cxxopts::value< uint32_t >()->default_value("90"))
			("m,match_value", "Smith-Waterman Match Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("n,sample_limit", "If the number of reads exceed this number then Graphite will randomly sample n reads [optional - default is no sample limit]", cxxopts::value< int32_t >()->default_value("-1"))
			("x,mismatch_value", "Smith-Waterman MisMatch Value [optional - default is 4]", cxxopts::value< uint32_t >()->default_value("4"))
			("g,gap_open_value", "Smith-Waterman Gap Open Value [optional - default is 6]", cxxopts::value< uint32_t >()->default_value("6"))
			("e,gap_extionsion_value", "Smith-Waterman Gap Extension Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("q,mapping_quality", "Mapping Quality Filter - (0 - 255) Filter reads that are less than or equal to this value [optional - default is no filter (-1)]", cxxopts::value< int32_t >()->default_value("-1"))
			("i,igv_visualization_output", "Output IGV input for visualization [optional - default is false]");
		this->m_options.parse(argc, argv);
	}

	void Params::parsePathTrace(int argc, char** argv)
	{

	}

	bool Params::showHelp()
	{
		return m_options.count("h");
	}

	void Params::printHelp()
	{
		std::cout << "Graphite Version: V" << GRAPHITE_VERSION_MAJOR << "." << GRAPHITE_VERSION_MINOR << "." << GRAPHITE_VERSION_PATCH << std::endl;
		std::cout << this->m_options.help() << std::endl;
	}

	bool Params::validateRequired()
	{
		std::vector< std::string > errorMessages;
		if (!m_options.count("v"))
		{
			errorMessages.emplace_back("vcf path(s) required");
		}
		if (!m_options.count("b"))
		{
			errorMessages.emplace_back("bam path(s) required");
		}
		if (!m_options.count("f"))
		{
			errorMessages.emplace_back("fasta path required");
		}
		if (!m_options.count("o"))
		{
			errorMessages.emplace_back("output path required");
		}
		if (m_options["q"].as< int32_t >() < -1 || m_options["q"].as< int32_t >() > 255)
		{
			errorMessages.emplace_back("invalid mapping quality, please proved a value between 0 and 255");
		}
		if (errorMessages.size() > 0)
		{
			std::cout << "There was a problem parsing commands" << std::endl;
			for (auto message : errorMessages)
			{
				std::cout << message << std::endl;
			}
			return false;
		}
		else
		{
			return true;
		}
	}


    bool Params::getIncludeDuplicates()
    {
        return m_options.count("d") > 0;
    }

	std::string Params::getFastaPath()
	{
		auto fastaPath = m_options["f"].as< std::string >();
		std::vector< std::string > fastaPaths = {fastaPath};
		validateFilePaths(fastaPaths, true);
		return fastaPath;
	}

	std::vector< std::string > Params::getInVCFPaths()
	{
		auto vcfPaths =  m_options["v"].as< std::vector< std::string > >();
		validateFilePaths(vcfPaths, true);
		return vcfPaths;
	}

	std::vector< std::string > Params::getBAMPaths()
	{
		auto bamPaths =  m_options["b"].as< std::vector< std::string > >();
		validateFilePaths(bamPaths, true);
		return bamPaths;
	}

	std::string Params::getOutputDirectory()
	{
		auto outputDir = m_options["o"].as< std::string >();
		std::vector< std::string > outputDirs = {outputDir};
		validateFolderPaths(outputDirs, true);
		return outputDir;
	}

	std::string Params::getOverwrittenSampleName()
	{
		return m_options["s"].as< std::string >();
	}

	Region::SharedPtr Params::getRegion()
	{
		if (m_options.count("r"))
		{
			return std::make_shared< Region >(m_options["r"].as< std::string >(), Region::BASED::ONE);
		}
		return nullptr;
	}

	uint32_t Params::getThreadCount()
	{
		return m_options["t"].as< uint32_t >();
	}

	int Params::getMatchValue()
	{
		return m_options["m"].as< uint32_t >();
	}

	int Params::getMisMatchValue()
	{
		return m_options["x"].as< uint32_t >();
	}

	int Params::getGapOpenValue()
	{
		return m_options["g"].as< uint32_t >();
	}

	int32_t Params::getMappingQualityFilter()
	{
		return m_options["q"].as< int32_t >();
	}

	int Params::getGapExtensionValue()
	{
		return m_options["e"].as< uint32_t >();
	}

	bool Params::outputVisualizationFiles()
	{
		return m_options["i"].as< bool >();
	}

	void Params::validateFolderPaths(const std::vector< std::string >& paths, bool exitOnFailure)
	{
		for (auto path : paths)
		{
			folderExists(path, exitOnFailure);
		}
	}

	void Params::validateFilePaths(const std::vector< std::string >& paths, bool exitOnFailure)
	{
		for (auto path : paths)
		{
			fileExists(path, exitOnFailure);
		}
	}
	int32_t Params::getReadSampleNumber()
	{
		return m_options["n"].as< int32_t >();
	}

}
