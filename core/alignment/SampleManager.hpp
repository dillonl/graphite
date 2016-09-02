#ifndef GRAPHITE_SAMPLE_MANAGER_HPP
#define GRAPHITE_SAMPLE_MANAGER_HPP

#include "Sample.hpp"
#include "core/util/Noncopyable.hpp"

#include <unordered_map>

namespace graphite
{
	class SampleManager : private Noncopyable
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
			std::lock_guard< std::mutex > lock(this->m_sample_ptrs_lock);
			auto iter = m_sample_ptrs_map.find(readGroup);
			if (iter != m_sample_ptrs_map.end())
			{
				return iter->second;
			}
			return nullptr;
		}

		void addSamplePtr(Sample::SharedPtr samplePtr)
		{
			std::lock_guard< std::mutex > lock(this->m_sample_ptrs_lock);
			if (m_sample_ptrs_map.find(samplePtr->getReadgroup()) == m_sample_ptrs_map.end())
			{
				m_sample_ptrs_map.emplace(samplePtr->getReadgroup(), samplePtr);
			}
		}

	private:
		SampleManager() {}
		~SampleManager() {}
		std::mutex m_sample_ptrs_lock;
		std::unordered_map< std::string, Sample::SharedPtr > m_sample_ptrs_map;
	};
}

#endif //GRAPHITE_SAMPLE_MANAGER_HPP
