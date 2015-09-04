#ifndef GRAPHITE_MAPPINGALIGNMENT_INFO_H
#define GRAPHITE_MAPPINGALIGNMENT_INFO_H

#include "core/alignment/IAlignment.h"

#include "externals/gssw/gssw.h"

#include <boost/noncopyable.hpp>

namespace graphite
{
	class IAdjudicator;
	class MappingAlignment : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< MappingAlignment > SharedPtr;
	MappingAlignment(IAlignment::SharedPtr alignmentPtr, uint32_t mappingOffset, uint32_t mappingOffsetLength, uint32_t mappingScore, std::vector< char >& cigarOperations, std::vector< uint32_t >& cigarOperationLengths) :
		m_alignment_ptr(alignmentPtr), m_mapping_offset(mappingOffset), m_mapping_offset_length(mappingOffsetLength), m_mapping_score(mappingScore), m_cigar_operations(cigarOperations) , m_cigar_operation_lengths(cigarOperationLengths), m_mapping_sequence(std::string(this->m_alignment_ptr->getSequence() + this->m_mapping_offset, this->m_mapping_offset_length))
		{
		}

		~MappingAlignment() {}

		IAlignment::SharedPtr getAlignmentPtr() { return this->m_alignment_ptr; }
		uint32_t getMappingScore(std::shared_ptr< IAdjudicator > adjudicatorPtr)
		{
			/*
			if (this->m_mapping_score < 0) { this->m_mapping_score = 0; }
			gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
			for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
			{
				for (int j = 0; j < nc->cigar->length; ++j)
				{
					switch (nc->cigar->elements[j].type)
					{
					case 'M':
						this->m_mapping_score += (adjudicatorPtr->getMatchValue() * nc->cigar->elements[j].length);
						break;
					case 'X':
						this->m_mapping_score -= (adjudicatorPtr->getMisMatchValue() * nc->cigar->elements[j].length);
						break;
					case 'I': // I and D are treated the same
					case 'D':
						this->m_mapping_score -= adjudicatorPtr->getGapOpenValue();
						this->m_mapping_score -= (adjudicatorPtr->getGapExtensionValue() * (nc->cigar->elements[j].length -1));
						break;
					default:
					}
					this->m_mapping_score = (this->m_mapping_score < 0) ? 0 : this->m_mapping_score; // the floor of the mapping score is 0
				}
			}
			*/
			return this->m_mapping_score;
		}
		uint32_t getMappingOffset() { return this->m_mapping_offset; }
		uint32_t getMappingOffsetLength() { return this->m_mapping_offset_length; }
		std::string getMappingSequence() { return this->m_mapping_sequence; }
		std::vector< char > getCigarOperations() { return this->m_cigar_operations; }
		std::vector< uint32_t > getCigarOperationLengths() { return this->m_cigar_operation_lengths; }

	private:
		IAlignment::SharedPtr m_alignment_ptr;
		int32_t m_mapping_score;
		uint32_t m_mapping_offset;
		uint32_t m_mapping_offset_length;
		std::string m_mapping_sequence;
		std::vector< char > m_cigar_operations;
		std::vector< uint32_t > m_cigar_operation_lengths;
	};
}

#endif //GRAPHITE_MAPPINGALIGNMENT_INFO_H
