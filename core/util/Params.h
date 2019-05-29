#ifndef GRAPHITE_PARAMS_H
#define GRAPHITE_PARAMS_H

#include "Noncopyable.hpp"
#include "core/region/Region.h"

#include <cstring>
#include <vector>

#include <cxxopts.hpp>


namespace graphite
{
	class Params : private Noncopyable
	{
	public:
		Params();
		~Params();

		void parseGSSW(int argc, char** argv);
		void parsePathTrace(int argc, char** argv);
		bool showHelp();
		void printHelp();
		bool validateRequired();

		std::string getFastaPath();
		std::vector< std::string > getInVCFPaths();
		std::vector< std::string > getBAMPaths();
		std::string getOutputDirectory();
		Region::SharedPtr getRegion();
        bool getIncludeDuplicates();
		uint32_t getPercent();
		uint32_t getThreadCount();
		int getMatchValue();
		int getMisMatchValue();
		int getGapOpenValue();
		int getGapExtensionValue();
		uint32_t getGraphSize();
		int32_t getMappingQualityFilter();
		bool outputVisualizationFiles();
		int32_t getReadSampleNumber();
	private:
		void validateFolderPaths(const std::vector< std::string >& paths, bool exitOnFailure);
		void validateFilePaths(const std::vector< std::string >& paths, bool exitOnFailure);

		cxxopts::Options m_options;

		std::vector< std::string > m_bam_paths;
	};
}

#endif //GRAPHITE_PARAMS_H
