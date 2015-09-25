#include "Params.h"

#include <thread>

namespace graphite
{
	Params::Params()
	{

	}

	Params::~Params()
	{
	}

	void Params::parseGSSW(int argc, char** argv)
	{
		m_options_description_ptr = std::make_shared< boost::program_options::options_description >("options");
		m_options_description_ptr->add_options()
			("help,h","Print help message")
			(",b", boost::program_options::value< std::string >()->required(), "Path to input BAM file")
			(",r", boost::program_options::value< std::string >()->required(), "Region information")
			(",v", boost::program_options::value< std::vector< std::string > >()->required()->multitoken(), "Path to input VCF file[s], separate multiple files by space")
			(",o", boost::program_options::value< std::string >()->default_value(""), "Path to output directory [optional - if not provided then prints to std::cout]")
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
	}

	void Params::parsePathTrace(int argc, char** argv)
	{
		m_options_description_ptr = std::make_shared< boost::program_options::options_description >("options");
		m_options_description_ptr->add_options()
			("help,h","Print help message")
			(",r", boost::program_options::value< std::string >()->required(), "Region information")
			(",v", boost::program_options::value< std::vector< std::string > >()->required()->multitoken(), "Path to input VCF file")
			(",f", boost::program_options::value< std::string >()->required(), "Path to input FASTA file")
		    (",p", boost::program_options::value< std::string >()->default_value(""), "Prefix to output files [optional - default is stdout]");
		auto parseCommandLine = boost::program_options::parse_command_line(argc, argv, *m_options_description_ptr);
		boost::program_options::store(parseCommandLine, m_variables_map);
	}

	bool Params::showHelp()
	{
		return (m_variables_map.count("help") > 0);
	}

	void Params::printHelp()
	{
		std::cout << *m_options_description_ptr << std::endl;
	}

	bool Params::validateRequired()
	{
		try
		{
			boost::program_options::notify(m_variables_map);
		}
		catch (boost::program_options::required_option& e)
		{
			std::cout << "ERROR " << e.what() << std::endl << std::endl;
			return false;
		}
		return true;
	}

	std::string Params::getFastaPath()
	{
		return m_variables_map["-f"].as< std::string >();
	}

	std::vector< std::string > Params::getInVCFPaths()
	{
		return m_variables_map["-v"].as< std::vector< std::string > >();
	}

	std::string Params::getBAMPath()
	{
		return m_variables_map["-b"].as< std::string >();
	}

	std::string Params::getOutputDirectory()
	{
		return m_variables_map["-o"].as< std::string >();
	}

	std::string Params::getFilePrefix()
	{
		return m_variables_map["-p"].as< std::string >();
	}

	Region::SharedPtr Params::getRegion()
	{
		auto regionString = m_variables_map["-r"].as< std::string >();
		return std::make_shared< Region >(regionString);
	}

	uint32_t Params::getPercent()
	{
		return m_variables_map["-p"].as< uint32_t >();
	}

	uint32_t Params::getThreadCount()
	{
		return m_variables_map["-t"].as< uint32_t >();
	}

	uint32_t Params::getGraphSize()
	{
		return m_variables_map["-g"].as< uint32_t >();
	}

	int Params::getMatchValue()
	{
		return m_variables_map["-m"].as< uint32_t >();
	}

	int Params::getMisMatchValue()
	{
		return m_variables_map["-s"].as< uint32_t >();
	}

	int Params::getGapOpenValue()
	{
		return m_variables_map["-a"].as< uint32_t >();
	}

	int Params::getGapExtensionValue()
	{
		return m_variables_map["-e"].as< uint32_t >();
	}

}
