#include "Sample.h"

namespace graphite
{
	Sample::Sample(const std::string& sampleName, const std::string& readGroup, const std::string& samplePath) :
		m_sample_name(sampleName),
		m_sample_readgroup(readGroup),
		m_sample_path(samplePath)
	{
	}

	Sample::~Sample()
	{
	}

	std::string Sample::getName() { return m_sample_name; }
	std::string Sample::getReadgroup() { return m_sample_readgroup; }
	std::string Sample::getPath() { return m_sample_path; }
}
