#ifndef GRAPHITE_IALIGNMENTREADER_HPP
#define GRAPHITE_IALIGNMENTREADER_HPP

#include "IAlignmentList.h"
#include "core/region/Region.h"
#include "core/util/Noncopyable.hpp"
#include "core/sample/SampleManager.h"

#include <memory>

namespace graphite
{
	class IAlignmentReader : private Noncopyable
	{
	public:
		IAlignmentReader()
		{
			static uint32_t s_id = 0;
			m_id = s_id++;
		}

		virtual ~IAlignmentReader() {}

		virtual void open() = 0;
		virtual void close() = 0;

		virtual std::vector< IAlignment::SharedPtr > loadAlignmentsInRegion(Region::SharedPtr regionPtr, SampleManager::SharedPtr sampleManagerPtr, bool unmappedOnly, bool excludeDuplicateReads) = 0;

		uint32_t getID() { return m_id; }

	protected:
		uint32_t m_id;
	};
}

#endif //GRAPHITE_IALIGNMENTREADER_HPP
