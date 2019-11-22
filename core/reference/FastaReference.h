#ifndef GRAPHITE_FASTAREFERENCE_H
#define GRAPHITE_FASTAREFERENCE_H

#include "core/util/Noncopyable.hpp"
#include "core/region/Region.h"

#include "Fasta.h"

#include <memory>

namespace graphite
{
	class FastaReference : private Noncopyable
	{
	public:
		typedef std::shared_ptr< FastaReference > SharedPtr;
        FastaReference(const std::string& fastaPath);
		~FastaReference();

		std::string getSequenceStringFromRegion(Region::SharedPtr regionPtr);

	private:
		void setReferenceIDAndSequence(Region::SharedPtr regionPtr);
		std::shared_ptr< ::FastaReference > m_fasta_reference;
		std::string m_current_reference_id;
		std::string m_reference_sequence;
	};
}

#endif //GRAPHITE_FASTAREFERENCE_H
