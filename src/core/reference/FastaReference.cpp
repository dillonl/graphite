#include "FastaReference.h"

namespace gwiz
{

	FastaReference::FastaReference(std::string path) :
		m_fasta_path(path)
	{
		m_fasta_reference = std::make_shared< fastahack::FastaReference >();
		m_fasta_reference->open(this->m_fasta_path);
	}

	FastaReference::~FastaReference()
	{
	}

	const char* FastaReference::getReferenceAtRegion(Region::SharedPtr region)
	{
		if (this->m_region_mapped_reference.find(region->getRegionString()) == this->m_region_mapped_reference.end())
		{
			std::string seqName = region->getReferenceID();
			position startPosition = region->getStartPosition();
			size_t length = region->getEndPosition() - startPosition;
			this->m_region_mapped_reference[region->getRegionString()] = this->m_fasta_reference->getSubSequence(seqName, startPosition, length);
		}
		return this->m_region_mapped_reference[region->getRegionString()].c_str();
	}

} // end namespace gwiz
