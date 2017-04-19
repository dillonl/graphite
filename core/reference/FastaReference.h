#ifndef GRAPHITE_FASTA_REFERENCE_H
#define GRAPHITE_FASTA_REFERENCE_H

#include "IReference.h"

#include "core/region/Region.h"
#include "Fasta.h"

#include <map>
#include <memory>

namespace graphite
{
	class FastaReference : public IReference
	{
	public:
		FastaReference(const std::string& path, Region::SharedPtr region);
		~FastaReference();

		const char* getSequence() override { return m_sequence.c_str(); }
		size_t getSequenceSize() override { return m_sequence.size();  }

	private:
		void setSequence(Region::SharedPtr region);

		std::string m_fasta_path;
		std::shared_ptr< ::FastaReference > m_fasta_reference;

	};
} // end namespace graphite

#endif // GRAPHITE_FASTA_REFERENCE_H
