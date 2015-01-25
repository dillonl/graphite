#ifndef GWIZ_FASTA_REFERENCE_H
#define GWIZ_FASTA_REFERENCE_H

#include "IReference.h"

#include "core/region/Region.h"
#include "fastahack/Fasta.h"

#include <map>
#include <memory>

namespace gwiz
{
	class FastaReference : public IReference
	{
	public:
		FastaReference(const std::string& path, Region::SharedPtr region);
		~FastaReference();

		const char* getSequence() { return m_sequence.c_str(); }

	private:
		void setRegion();

		std::string m_fasta_path;
		Region::SharedPtr m_region;
		std::string m_sequence;
		std::shared_ptr< fastahack::FastaReference > m_fasta_reference;

	};
} // end namespace gwiz

#endif // GWIZ_FASTA_REFERENCE_H
