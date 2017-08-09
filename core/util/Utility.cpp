#include "Utility.h"

#include "core/file/ASCIIFileWriter.h"
#include "core/file/BGZFFileWriter.h"

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

	std::unordered_map< std::string, IFileWriter::SharedPtr > getUniqueFileNames(const std::vector< std::string >& filePaths, const std::string& outputDirectory)
	{
		std::unordered_map< std::string, IFileWriter::SharedPtr > outPaths;
		for (auto filePath : filePaths)
		{
			std::string path = filePath.substr(filePath.find_last_of("/") + 1);
			std::string outFilePath = outputDirectory + "/" + path;
			uint32_t counter = 1;
			std::string extension = "";
			while (IFile::fileExists(outFilePath, false))
			{
				extension = filePath.substr(filePath.find_last_of(".") + 1);
				std::string fileNameWithoutExtension = path.substr(0, path.find_last_of("."));
				outFilePath = outputDirectory + "/" + fileNameWithoutExtension + "." + std::to_string(counter) + "." + extension;
				++counter;
			}
			IFileWriter::SharedPtr fileWriterPtr;
			if (extension.compare("gz") == 0)
			{
				fileWriterPtr = std::make_shared< BGZFFileWriter >(outFilePath);
			}
			else
			{
				fileWriterPtr = std::make_shared< ASCIIFileWriter >(outFilePath);
			}
			fileWriterPtr->open();
			outPaths[filePath] = fileWriterPtr;
		}
		return outPaths;
	}
}
