#ifndef GWIZ_FASTA_REFERENCE_H
#define GWIZ_FASTA_REFERENCE_H

#include "IReference.h"

#include "fastahack/Fasta.h"

#include <map>
#include <memory>

namespace gwiz
{
	class FastaReference : IReference
	{
	public:
		FastaReference(std::string path);
		~FastaReference();

		const char* getReferenceAtRegion(Region::SharedPtr region) override;
	private:
		std::string m_fasta_path;
		std::map< std::string, std::string > m_region_mapped_reference; // key = region string, value = reference at that region

		std::shared_ptr< fastahack::FastaReference > m_fasta_reference;

	};
} // end namespace gwiz

#endif // GWIZ_FASTA_REFERENCE_H
