#ifndef GWIZ_GENOTYPERALLELE_H
#define GWIZ_GENOTYPERALLELE_H

#include "core/alignment/IAlignment.h"

#include <memory>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class GenotyperAllele : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< GenotyperAllele > SharedPtr;

		enum class Type {REFERENCE, VARIANT_ALTERNATE, VARIANT_REFERENCE};

		static std::string TypeToString(Type type)
		{
			switch (type)
			{
			case Type::VARIANT_ALTERNATE:
				return "VA";
			case Type::VARIANT_REFERENCE:
				return "VR";
			default:
				return "R";
			}
		}

	    GenotyperAllele(const Type alleleType, const std::string& sequence, position pos) :
		    m_read_count(0),
			m_allele_type(alleleType),
			m_position(pos),
			m_sequence(sequence)
		{
		}

		virtual ~GenotyperAllele()
		{
		}

		std::string getSequence() const { return m_sequence; }
		// void incrementReadCount() { ++this->m_read_count; }
		// uint32_t getReadCount() { return this->m_read_count; }
		uint32_t getReadCount() { return this->m_alignments.size(); }
		position getPosition() {return this->m_position; }
		Type getType() { return this->m_allele_type; }
		void addAlignment(IAlignment::SharedPtr alignmentPtr) { this->m_alignments.push_back(alignmentPtr); }
		std::list< IAlignment::SharedPtr > getAlignments() { return this->m_alignments; }

	protected:
		std::list< IAlignment::SharedPtr > m_alignments;
		uint32_t m_read_count;
		Type m_allele_type;
		std::string m_sequence;
		position m_position;
	};
}

#endif //GWIZ_GENOTYPERALLELE_H
