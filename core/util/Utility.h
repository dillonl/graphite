#ifndef GRAPHITE_CORE_UTIL_UTILITY_H
#define GRAPHITE_CORE_UTIL_UTILITY_H

#include "core/file/IFile.h"
#include "core/file/IFileWriter.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace graphite
{
	void split(const std::string& s, char c, std::vector< std::string >& v);

	std::unordered_map< std::string, IFileWriter::SharedPtr > getUniqueFileNames(const std::vector< std::string >& filePaths, const std::string& outputDirectory);
}

#endif //GRAPHITE_CORE_UTIL_UTILITY_H
