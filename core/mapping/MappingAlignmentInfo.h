#ifndef GRAPHITE_MAPPINGALIGNMENT_INFO_H
#define GRAPHITE_MAPPINGALIGNMENT_INFO_H

#include "core/alignment/IAlignment.h"

//#include "gssw.h"
#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class IAdjudicator;
	class MappingAlignmentInfo : private Noncopyable
	{
	public:
		typedef std::shared_ptr< MappingAlignmentInfo > SharedPtr;
	    MappingAlignmentInfo(IAllele::SharedPtr allelePtr, uint32_t swScore, uint32_t length, uint32_t prefixMatch, uint32_t suffixMatch) :
		    m_allele_ptr(allelePtr), m_sw_score(swScore), m_length(length), m_prefix_match(prefixMatch), m_suffix_match(suffixMatch)
		{
		}
		~MappingAlignmentInfo() {}

		uint32_t getSWScore() { return m_sw_score; }
		uint32_t getLength() { return m_length; }
		uint32_t getPrefixMatch() { return m_prefix_match; }
		uint32_t getSuffixMatch() { return m_suffix_match; }
		IAllele::SharedPtr getAllelePtr() { return m_allele_ptr; }
	private:
		uint32_t m_sw_score;
		uint32_t m_length;
		uint32_t m_prefix_match;
		uint32_t m_suffix_match;
		IAllele::SharedPtr m_allele_ptr;
	};
}

#endif //GRAPHITE_MAPPINGALIGNMENT_INFO_H
