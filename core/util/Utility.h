#ifndef GRAPHITE_CORE_UTIL_UTILITY_H
#define GRAPHITE_CORE_UTIL_UTILITY_H

#include <string>
#include <vector>

namespace graphite
{
	void split(const std::string& s, char c, std::vector< std::string >& v);
	bool fileExists(const std::string& name, bool exitOnFailure);
	bool folderExists(const std::string& path, bool exitOnFailure);
}

#endif //GRAPHITE_CORE_UTIL_UTILITY_H
