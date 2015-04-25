#ifndef GWIZ_PARAMETERS_H
#define GWIZ_PARAMETERS_H

#include <boost/noncopyable.hpp>

#include <string>
#include <map>

namespace gwiz
{

	class Parameters : private boost::noncopyable
	{
	public:
		static Parameters* Instance();

		std::string getParameterValue(const std::string& key);
		void setParams(const int argc, char** argv);
		void printUsage();

		std::string getBAMPath() { return m_bam_path; }
		std::string getFastaPath() { return m_fasta_path; }
		std::string getInVCFPath() { return m_in_vcf_path; }
		std::string getOutVCFPath() { return m_out_vcf_path; }
		std::string getRegion() { return m_region; }

	private:
		Parameters();
		~Parameters();
		bool setPath(std::string& outPath, const std::string& inPath);

		std::map<std::string, std::string> m_params;
		std::string m_command_options;

		std::string m_bam_path;
		std::string m_fasta_path;
		std::string m_in_vcf_path;
		std::string m_out_vcf_path;
		std::string m_region;

		static Parameters* s_instance;
	};

} // end namespace gwiz

#endif // BAMQ_PARAMETERS_H
