#ifndef GRAPHITE_SAMPLE_MANAGER_HPP
#define GRAPHITE_SAMPLE_MANAGER_HPP

#include "Sample.hpp"

#include <boost/noncopyable.hpp>
#include <unordered_map>

namespace graphite
{
	class SampleManager : private boost::noncopyable
	{
	public:
		static SampleManager* Instance()
		{
			static SampleManager* s_instance = nullptr; // lazy loaded
			if (s_instance == nullptr)
			{
				s_instance = new SampleManager();
			}
			return s_instance;
		}

		Sample::SharedPtr getSamplePtr(const std::string& readGroup)
		{
			auto iter = m_sample_ptrs_map.find(readGroup);
			// std::cout << "fd: " << *iter << std::endl;
			if (iter != m_sample_ptrs_map.end())
			{
				return iter->second;
			}
			return nullptr;
		}

		void addSamplePtr(Sample::SharedPtr samplePtr)
		{
			if (m_sample_ptrs_map.find(samplePtr->getReadgroup()) == m_sample_ptrs_map.end())
			{
				m_sample_ptrs_map.emplace(samplePtr->getReadgroup(), samplePtr);
			}
		}

	private:
		std::unordered_map< std::string, Sample::SharedPtr > m_sample_ptrs_map;
	};
}

#endif //GRAPHITE_SAMPLE_MANAGER_HPP
