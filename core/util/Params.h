#ifndef GWIZ_PARAMS_H
#define GWIZ_PARAMS_H

#include "core/region/Region.h"
#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>

namespace gwiz
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
		std::string getBAMPath();
		std::string getOutVCFPath();
		std::string getFilePrefix();
		Region::SharedPtr getRegion();
		uint32_t getPercent();
		uint32_t getThreadCount();
		int getMatchValue();
		int getMisMatchValue();
		int getGapOpenValue();
		int getGapExtensionValue();
	private:
		std::shared_ptr< boost::program_options::options_description > m_options_description_ptr;
		boost::program_options::variables_map m_variables_map;
	};
}

#endif //GWIZ_PARAMS_H
