#include "Params.h"
#include "Utility.h"

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
			("r,region", "Region information", cxxopts::value< std::string >()->default_value(""))
			("o,output_directory", "Path to output directory", cxxopts::value< std::string >())
			("f,fasta", "Path to input FASTA file", cxxopts::value< std::string >())
			("p,percent_match", "Smith-Waterman Percent [optional - default is 90]", cxxopts::value< uint32_t >()->default_value("90"))
			("m,match_value", "Smith-Waterman Match Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("s,mismatch_value", "Smith-Waterman MisMatch Value [optional - default is 4]", cxxopts::value< uint32_t >()->default_value("4"))
			("a,gap_open_value", "Smith-Waterman Gap Open Value [optional - default is 6]", cxxopts::value< uint32_t >()->default_value("6"))
			("e,gap_extionsion_value", "Smith-Waterman Gap Extension Value [optional - default is 1]", cxxopts::value< uint32_t >()->default_value("1"))
			("g,graph_size", "The size of the graph [optional - default is 3000]", cxxopts::value< uint32_t >()->default_value("3000"))
			("t,number_threads", "Thread count [optional - default is number of cores x 2]", cxxopts::value< uint32_t >()->default_value(std::to_string(std::thread::hardware_concurrency() * 2)));
		this->m_options.parse(argc, argv);
		/*
		m_options_description_ptr = std::make_shared< boost::program_options::options_description >("options");
		m_options_description_ptr->add_options()
			("help,h","Print help message")
			(",b", boost::program_options::value< std::vector< std::string > >()->required()->multitoken(), "Path to input BAM file[s], separate multiple files by space")
			(",r", boost::program_options::value< std::string >()->default_value(""), "Region information")
			(",d", boost::program_options::bool_switch()->default_value(false), "Exclude Duplicate Reads")
			(",v", boost::program_options::value< std::vector< std::string > >()->required()->multitoken(), "Path to input VCF file[s], separate multiple files by space")
			(",o", boost::program_options::value< std::string >()->required(), "Path to output directory")
			(",f", boost::program_options::value< std::string >()->required(), "Path to input FASTA file")
			(",p", boost::program_options::value< uint32_t >()->default_value(90), "Smith-Waterman Percent [optional - default is 90]")
			(",m", boost::program_options::value< uint32_t >()->default_value(1), "Smith-Waterman Match Value [optional - default is 1]")
			(",s", boost::program_options::value< uint32_t >()->default_value(4), "Smith-Waterman MisMatch Value [optional - default is 4]")
			(",a", boost::program_options::value< uint32_t >()->default_value(6), "Smith-Waterman Gap Open Value [optional - default is 6]")
			(",e", boost::program_options::value< uint32_t >()->default_value(1), "Smith-Waterman Gap Extension Value [optional - default is 1]")
			(",g", boost::program_options::value< uint32_t >()->default_value(3000), "The size of the graph [optional - default is 3000]")
			(",t", boost::program_options::value< uint32_t >()->default_value(std::thread::hardware_concurrency() * 2), "Thread count [optional - default is number of cores x 2]");
		auto parseCommandLine = boost::program_options::parse_command_line(argc, argv, *m_options_description_ptr);
		boost::program_options::store(parseCommandLine, m_variables_map);
		*/
	}

	void Params::parsePathTrace(int argc, char** argv)
	{
		/*
		m_options_description_ptr = std::make_shared< boost::program_options::options_description >("options");
		m_options_description_ptr->add_options()
			("help,h","Print help message")
			(",r", boost::program_options::value< std::string >()->required(), "Region information")
			(",v", boost::program_options::value< std::vector< std::string > >()->required()->multitoken(), "Path to input VCF file")
			(",f", boost::program_options::value< std::string >()->required(), "Path to input FASTA file")
		    (",p", boost::program_options::value< std::string >()->default_value(""), "Prefix to output files [optional - default is stdout]");
		auto parseCommandLine = boost::program_options::parse_command_line(argc, argv, *m_options_description_ptr);
		boost::program_options::store(parseCommandLine, m_variables_map);
		*/
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
		return m_options["f"].as< std::string >();
	}

	std::vector< std::string > Params::getInVCFPaths()
	{
		return m_options["v"].as< std::vector< std::string > >();
	}

	std::vector< std::string > Params::getBAMPaths()
	{
		return m_options["b"].as< std::vector< std::string > >();
	}

	std::string Params::getOutputDirectory()
	{
		return m_options["o"].as< std::string >();
	}

	std::string Params::getFilePrefix()
	{
		// return m_variables_map["-p"].as< std::string >();
		return "";
	}

	Region::SharedPtr Params::getRegion()
	{

		if (m_options.count("r"))
		{
			return std::make_shared< Region >(m_options["r"].as< std::string >());
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

}
