#include "Parameters.h"

#include <getopt.h>
#include <iostream>

namespace gwiz
{
	Parameters* Parameters::s_instance = NULL;

	Parameters::Parameters() :
		m_command_options("h:b:f:v:o:r:"),
		m_bam_path(""),
		m_in_vcf_path(""),
		m_out_vcf_path(""),
		m_region("")
	{
	}

	Parameters::~Parameters()
	{
	}

	Parameters* Parameters::Instance()
	{

		if (s_instance == NULL)
		{
			s_instance = new Parameters();
		}
		return s_instance;
	}

	std::string Parameters::getParameterValue(const std::string& key)
	{
		if (m_params.find(key) == m_params.end())
		{
			return "";
		}
		return m_params[key];
	}

	void Parameters::setParams(const int argc, char** argv)
	{
		m_params.clear();

		int option_index = 0;
		static struct option long_options[] =
			{
				{"help", no_argument, 0, 'h'},
				{"bam", required_argument, 0, 'b'},
				{"fasta", required_argument, 0, 'f'},
				{"vcf", required_argument, 0, 'v'},
				{"output_vcf", required_argument, 0, 'o'},
				{"region", required_argument, 0, 'r'},
				{NULL, 0, 0, 0}
			};
		char c;
		while ((c = getopt_long(argc, argv, m_command_options.c_str(), long_options, NULL)) != -1)
		{
			switch (c)
			{
			case 'h':
				printUsage();
				exit(0);
			case 'f':
				if (!setPath(this->m_fasta_path, optarg))
				{
					std::cout << "Invalid Fasta Path: " << optarg << std::endl;
					exit(0);
				}
				break;
			case 'b':
				if (!setPath(this->m_bam_path, optarg))
				{
					std::cout << "Invalid Bam Path: " << optarg << std::endl;
					exit(0);
				}
				break;
			case 'v':
				if (!setPath(this->m_in_vcf_path, optarg))
				{
					std::cout << "Invalid Input VCF Path: " << optarg << std::endl;
					exit(0);
				}
				break;
			case 'o':
				this->m_out_vcf_path = optarg;
				break;
			case 'r':
				this->m_region = optarg;
				break;
			default:
				std::cout << "Invalid Parameter" << std::endl;
				printUsage();
				exit(0);

			}
		}

	}

	void Parameters::printUsage()
	{
		std::cout << "Options:" << std::endl;
		std::cout << "\th\tPrints this statement." << std::endl;
		std::cout << "\tb\tPath to input BAM file" << std::endl;
		std::cout << "\tv\tPath to input VCF file" << std::endl;
		std::cout << "\to\tPath to output VCF file [optional - default is stdout]" << std::endl;
	}

	// A non-member function to check if a file exists
	inline bool DoesFileExist(const std::string& path)
	{
		if (FILE* file = fopen(path.c_str(), "r"))
		{
			fclose(file);
			return true;
		}
		else
		{
			return false;
		}
	}

	bool Parameters::setPath(std::string& outPath, const std::string& inPath)
	{
		if (DoesFileExist(inPath))
		{
			outPath = inPath;
			return true;
	    }
		else
		{
			return false;
		}
	}


} // end namespace gwiz
