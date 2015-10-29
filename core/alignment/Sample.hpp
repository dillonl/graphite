#ifndef GRAPHITE_SAMPLE_HPP
#define GRAPHITE_SAMPLE_HPP

#include <string>
#include <memory>

namespace graphite
{
	class Sample : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< Sample > SharedPtr;
		Sample(const std::string& sampleName, const std::string& readGroup, const std::string& samplePath) :
			m_sample_name(sampleName),
			m_sample_readgroup(readGroup),
			m_sample_path(samplePath)
		{
		}

		std::string getName() { return m_sample_name; }
		std::string getReadgroup() { return m_sample_readgroup; }
		std::string getPath() { return m_sample_path; }

	private:
		std::string m_sample_name;
		std::string m_sample_readgroup;
		std::string m_sample_path;
	};
}

#endif //GRAPHITE_SAMPLE_HPP
