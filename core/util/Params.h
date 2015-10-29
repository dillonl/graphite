#ifndef GRAPHITE_PARAMS_H
#define GRAPHITE_PARAMS_H

#include "core/region/Region.h"
#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>

namespace graphite
{
	class Params : private boost::noncopyable
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
		uint32_t getPercent();
		uint32_t getThreadCount();
		int getMatchValue();
		int getMisMatchValue();
		int getGapOpenValue();
		int getGapExtensionValue();
		uint32_t getGraphSize();
	private:
		std::shared_ptr< boost::program_options::options_description > m_options_description_ptr;
		boost::program_options::variables_map m_variables_map;
	};
}

#endif //GRAPHITE_PARAMS_H
