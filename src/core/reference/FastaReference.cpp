#include "FastaReference.h"

namespace gwiz
{

	FastaReference::FastaReference(const std::string& path, Region::SharedPtr region) :
		m_fasta_path(path), m_region(region)
	{
		m_fasta_reference = std::make_shared< fastahack::FastaReference >();
		m_fasta_reference->open(this->m_fasta_path);
		setRegion();
	}

	FastaReference::~FastaReference()
	{
	}

	void FastaReference::setRegion()
	{
		std::string seqName = this->m_region->getReferenceID();
		position startPosition = this->m_region->getStartPosition();
		size_t length = this->m_region->getEndPosition() - startPosition;
		this->m_sequence = this->m_fasta_reference->getSubSequence(seqName, startPosition, length);
	}

} // end namespace gwiz
