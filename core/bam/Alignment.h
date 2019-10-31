#pragma once

#include "core/util/Noncopyable.hpp"
#include "core/sample/Sample.h"

#include <string>
#include <memory>

namespace graphite
{
	class Alignment : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Alignment > SharedPtr;
		Alignment(char* htslibSeq, uint32_t len, std::string& readName, bool isForwardStrand, bool isFirstMate, uint16_t mapQuality, Sample::SharedPtr samplePtr);
		~Alignment();

		char* getSequence() { return this->m_seq; }
		uint32_t getLength() { return this->m_len; }
		std::string getReadName() { return this->m_read_name; }
		std::string getUniqueReadName() { return this->m_unique_read_name; }
		bool getIsForwardStrand() { return m_is_forward_strand; }
		bool getIsFirstMate() { return m_is_first_mate; }
		uint16_t getMapQuality() { return m_map_quality; }
		Sample::SharedPtr getSample() { return m_sample_ptr; }

	private:
		char* m_seq;
		uint32_t m_len;
		std::string m_read_name;
		std::string m_unique_read_name;
		bool m_is_forward_strand;
		bool m_is_first_mate;
		uint16_t m_map_quality;
		Sample::SharedPtr m_sample_ptr;
	};
}
