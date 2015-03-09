#ifndef GWIZ_GENOTYPERALLELE_H
#define GWIZ_GENOTYPERALLELE_H

#include <memory>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IGenotyperAllele : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IGenotyperAllele > SharedPtr;

		enum class Type {REFERENCE, ALTERNATE};

	    IGenotyperAllele(const Type alleleType, const std::string& sequence, const position pos) :
		    m_read_count(0),
			m_allele_type(alleleType),
			m_position(pos),
			m_sequence(sequence)
		{
		}

		virtual ~IGenotyperAllele()
		{
		}

		std::string getSequence() const { return m_sequence; }
		position getPosition() const { return m_position; }
		void incrementReadCount() { ++this->m_read_count; }
		uint32_t getReadCount() { return this->m_read_count; }


	protected:
		uint32_t m_read_count;
		Type m_allele_type;
		position m_position;
		std::string m_sequence;
	};
}

#endif //GWIZ_GENOTYPERALLELE_H
