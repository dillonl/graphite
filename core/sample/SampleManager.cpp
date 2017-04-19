#include "SampleManager.h"

namespace graphite
{
	SampleManager::SampleManager(const std::vector< Sample::SharedPtr >& samplePtrs)
	{
		for (auto samplePtr : samplePtrs)
		{
			addSamplePtr(samplePtr);
		}
	}

	SampleManager::~SampleManager() {}

	Sample::SharedPtr SampleManager::getSamplePtr(const std::string& readGroup)
	{
		std::lock_guard< std::mutex > lock(this->m_sample_ptrs_lock);
		auto iter = m_sample_ptrs_map.find(readGroup);
		if (iter != m_sample_ptrs_map.end())
		{
			return iter->second;
		}
		return nullptr;
	}

	void SampleManager::addSamplePtr(Sample::SharedPtr samplePtr)
	{
		std::lock_guard< std::mutex > lock(this->m_sample_ptrs_lock);
		if (m_sample_ptrs_map.find(samplePtr->getReadgroup()) == m_sample_ptrs_map.end())
		{
			m_sample_ptrs_map.emplace(samplePtr->getReadgroup(), samplePtr);
		}
	}

	uint32_t SampleManager::getSampleCount()
	{
		return m_sample_ptrs_map.size();
	}

	std::vector< Sample::SharedPtr > SampleManager::getSamplePtrs()
	{
		std::vector< Sample::SharedPtr > samplePtrs;
		for (auto iter : m_sample_ptrs_map)
		{
			samplePtrs.emplace_back(iter.second);
		}
		return samplePtrs;
	}
}
