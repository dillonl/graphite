#ifndef GWIZ_PARAMETERS_H
#define GWIZ_PARAMETERS_H

#include <string>
#include <map>

#include "utils/NonCopyable.h"

namespace gwiz
{

	class Parameters : private noncopyable
	{
	public:
		static Parameters* Instance();

		std::string getParameterValue(const std::string& key);
		void setParams(const int argc, char** argv);
		void printUsage();

		std::string getBamSourcePath() { return m_bam_source_path; }
		std::string getBamGliaPath() { return m_bam_glia_path; }
		std::string getRegion() { return m_region; }

	private:
		Parameters();
		~Parameters();

		bool setBamSourcePath(const std::string& path);
		bool setBamGliaPath(const std::string& path);

		std::map<std::string, std::string> m_params;
		std::string m_command_options;

		std::string m_bam_source_path;
		std::string m_bam_glia_path;
		std::string m_region;

		static Parameters* s_instance;
	};

} // end namespace gwiz

#endif // BAMQ_PARAMETERS_H
