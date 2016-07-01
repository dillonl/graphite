#ifndef GRAPHITE_SAMPLE_HPP
#define GRAPHITE_SAMPLE_HPP

#include <iostream>
#include <string>
#include <memory>

#include <boost/noncopyable.hpp>

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
		std::string getFoo() { return ""; }

	private:
		std::string m_sample_name;
		std::string m_sample_readgroup;
		std::string m_sample_path;
	};
}

#endif //GRAPHITE_SAMPLE_HPP
