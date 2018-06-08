#ifndef GRAPHITE_SAMPLE_MANAGER_H
#define GRAPHITE_SAMPLE_MANAGER_H

#include "Sample.h"
#include "core/util/Noncopyable.hpp"

#include <unordered_map>
#include <mutex>
#include <memory>
#include <vector>

namespace graphite
{
	class SampleManager : private Noncopyable
	{
	public:
		typedef std::shared_ptr< SampleManager > SharedPtr;
		SampleManager(const std::vector< Sample::SharedPtr >& samplePtrs);
		SampleManager(const std::vector< std::string >& bamPaths);
		~SampleManager();

		Sample::SharedPtr getSamplePtr(const std::string& readGroup);
		void addSamplePtr(Sample::SharedPtr samplePtr);
		uint32_t getSampleCount();
		std::vector< Sample::SharedPtr > getSamplePtrs();

	private:
		std::mutex m_sample_ptrs_lock;
		std::unordered_map< std::string, Sample::SharedPtr > m_sample_ptrs_map;
	};
}

#endif //GRAPHITE_SAMPLE_MANAGER_H
