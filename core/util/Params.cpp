#include "Params.h"
#include "Utility.h"
#include "core/file/IFile.h"

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
			("d,exclude_duplicates", "Exclude Duplicate Reads")
			("v,vcf", "Path to input VCF file[s], separate multiple files by space", cxxopts::value< std::vector< std::string > >())
			("b,bam", "Path to input BAM file[s], separate multiple files by space", cxxopts::value< std::vector< std::string > >())
			("r,region", "Region information", cxxopts::value< std::string >())
			("o,output_directory", "Path to output directory", cxxopts::value< std::string >())
			("f,fasta", "Path to input FASTA file", cxxopts::value< std::string >())
			("p,percent_match", "Smith-Waterman Percent [optional - default is 90]", cxxopts::value< uint32_t >()->default_value("90"))
			("m,match_value", "Smith-Waterman Match Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("s,mismatch_value", "Smith-Waterman MisMatch Value [optional - default is 4]", cxxopts::value< uint32_t >()->default_value("4"))
			("a,gap_open_value", "Smith-Waterman Gap Open Value [optional - default is 6]", cxxopts::value< uint32_t >()->default_value("6"))
			("e,gap_extionsion_value", "Smith-Waterman Gap Extension Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("g,graph_size", "The size of the graph [optional - default is 3000]", cxxopts::value< uint32_t >()->default_value("3000"))
			("t,number_threads", "Thread count [optional - default is number of cores x 2]", cxxopts::value< uint32_t >()->default_value(std::to_string(std::thread::hardware_concurrency() * 2)))
			("i,igv_files", "Output new fasta, bam and bed files for viewing in IGV");
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


    bool Params::getExcludeDuplicates()
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

	Region::SharedPtr Params::getRegion()
	{
		if (m_options.count("r"))
		{
			return std::make_shared< Region >(m_options["r"].as< std::string >(), Region::BASED::ONE);
		}
		return nullptr;
	}

	uint32_t Params::getPercent()
	{
		return m_options["p"].as< uint32_t >();
	}

	uint32_t Params::getThreadCount()
	{
		return m_options["t"].as< uint32_t >();
	}

	uint32_t Params::getGraphSize()
	{
		return m_options["g"].as< uint32_t >();
	}

	int Params::getMatchValue()
	{
		return m_options["m"].as< uint32_t >();
	}

	int Params::getMisMatchValue()
	{
		return m_options["s"].as< uint32_t >();
	}

	int Params::getGapOpenValue()
	{
		return m_options["a"].as< uint32_t >();
	}

	int Params::getGapExtensionValue()
	{
		return m_options["e"].as< uint32_t >();
	}

	bool Params::getIGVOutput()
	{
        return m_options.count("i") > 0;
	}

	void Params::validateFolderPaths(const std::vector< std::string >& paths, bool exitOnFailure)
	{
		for (auto path : paths)
		{
			IFile::folderExists(path, exitOnFailure);
		}
	}

	void Params::validateFilePaths(const std::vector< std::string >& paths, bool exitOnFailure)
	{
		for (auto path : paths)
		{
			IFile::fileExists(path, exitOnFailure);
		}
	}

}
