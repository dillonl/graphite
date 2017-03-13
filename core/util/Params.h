#ifndef GRAPHITE_PARAMS_H
#define GRAPHITE_PARAMS_H

#include "Noncopyable.hpp"
#include "core/region/Region.h"
#include <cxxopts.hpp>

#include <stdio.h>
#include <string.h>
#include <cstring>
#include <vector>

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
		std::string getFilePrefix();
		Region::SharedPtr getRegion();
        bool getExcludeDuplicates();
		uint32_t getPercent();
		uint32_t getThreadCount();
		int getMatchValue();
		int getMisMatchValue();
		int getGapOpenValue();
		int getGapExtensionValue();
		uint32_t getGraphSize();
	private:
		cxxopts::Options m_options;

		std::vector< std::string > m_bam_paths;
	};
}

#endif //GRAPHITE_PARAMS_H
