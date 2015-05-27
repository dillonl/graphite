#ifndef GWIZ_IALLELE_H
#define GWIZ_IALLELE_H

#include <core/sequence/Sequence.h>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class IAllele : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< IAllele > SharedPtr;
	IAllele(Sequence::SharedPtr sequencePtr) :
		m_sequence_ptr(sequencePtr) { }

		virtual ~IAllele() {}

		Sequence::SharedPtr getSequence() { return this->m_sequence_ptr; }

	protected:
		Sequence::SharedPtr m_sequence_ptr;

	private:

	};
}

#endif //GWIZ_IALLELE_H
