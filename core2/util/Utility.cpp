#include "Utility.h"

// #include <regex>
#include <string>
#include <iostream>

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


	/*
	void split(const std::string& s, std::vector< std::string >& v)
	{
		std::regex re("\\s+");
		std::sregex_token_iterator it(s.begin(), s.end(), re, -1);
		std::sregex_token_iterator reg_end;
		for (; it != reg_end; ++it)
		{
			std::cout << "s: " << it->str() << std::endl;
			v.emplace_back(it->str());
		}
	}
	*/
}
