#include "Utility.h"

// #include <regex>
#include <string>
#include <iostream>
#include <memory>
#include <fstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

namespace graphite
{
	void split(const std::string& s, char c, std::vector< std::string >& v)
	{
		std::string::size_type i = 0;
		std::string::size_type j = s.find(c);

		while (j != std::string::npos)
		{
			v.push_back(s.substr(i, j-i));
			i = ++j;
			j = s.find(c, j);

			if (j == std::string::npos)
			{
				v.push_back(s.substr(i, s.length()));
			}
		}
	}

	bool fileExists(const std::string& name, bool exitOnFailure)
	{
	    std::ifstream f(name.c_str());
	    if (!f)
	    {
			if (exitOnFailure)
			{
				std::cout << "File not found: " << name << std::endl;
				exit(EXIT_FAILURE);
			}
			return false;
		}
		return true;
	}

	bool folderExists(const std::string& path, bool exitOnFailure)
	{
		if(!path.empty())
		{
			struct stat sb;

			if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
			{
				return true;
			}
		}
		if (exitOnFailure)
		{
			std::cout << "Folder not found: " << path << std::endl;
			exit(EXIT_FAILURE);
		}
		return false;
	}

}
