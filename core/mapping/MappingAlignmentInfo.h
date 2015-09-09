#ifndef GRAPHITE_MAPPINGALIGNMENT_INFO_H
#define GRAPHITE_MAPPINGALIGNMENT_INFO_H

#include "core/alignment/IAlignment.h"

#include "externals/gssw/gssw.h"

#include <boost/noncopyable.hpp>

namespace graphite
{
	class IAdjudicator;
	class MappingAlignmentInfo : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< MappingAlignmentInfo > SharedPtr;
		MappingAlignmentInfo()
		{
		}

		void setSWScore(uint32_t swScore) { m_sw_score = swScore; }
		void setLength(uint32_t length) { m_length = length; }
		void setPrefixMatch(uint32_t prefixMatch) { m_prefix_match = prefixMatch; }
		void setSuffixMatch(uint32_t suffixMatch) { m_suffix_match = suffixMatch; }

		uint32_t getSWScore() { return m_sw_score; }
		uint32_t getLength() { return m_length; }
		uint32_t getPrefixMatch() { return m_prefix_match; }
		uint32_t getSuffixMatch() { return m_suffix_match; }

		~MappingAlignmentInfo() {}
	private:
		uint32_t m_sw_score;
		uint32_t m_length;
		uint32_t m_prefix_match;
		uint32_t m_suffix_match;
	};
}

#endif //GRAPHITE_MAPPINGALIGNMENT_INFO_H
