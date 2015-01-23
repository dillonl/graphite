#include "Parameters.h"

#include <getopt.h>
#include <iostream>

namespace gwiz
{
	Parameters* Parameters::s_instance = NULL;

	Parameters::Parameters() :
		m_command_options("hs:g:r:"),
		m_bam_source_path(""),
		m_bam_glia_path(""),
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
				{"bam_source", required_argument, 0, 's'},
				{"bam_glia", required_argument, 0, 'g'},
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
			case 's':
				if (!setBamSourcePath(optarg))
				{
					std::cout << "Invalid Bam Source Path: " << optarg << std::endl;
					exit(0);
				}
				break;
			case 'g':
				if (!setBamGliaPath(optarg))
				{
					std::cout << "Invalid Bam Glia Path: " << optarg << std::endl;
					exit(0);
				}
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

	bool Parameters::setBamSourcePath(const std::string& path)
	{
		if (DoesFileExist(path))
		{
			this->m_bam_source_path = path;
			return true;
	    }
		else
		{
			return false;
		}
	}

	bool Parameters::setBamGliaPath(const std::string& path)
	{
		if (DoesFileExist(path))
		{
			this->m_bam_glia_path = path;
			return true;
	    }
		else
		{
			return false;
		}
	}


} // end namespace gwiz
