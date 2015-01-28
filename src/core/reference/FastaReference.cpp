#include "FastaReference.h"

namespace gwiz
{

	FastaReference::FastaReference(const std::string& path, Region::SharedPtr region) :
		m_fasta_path(path),
		IReference()
	{
		m_fasta_reference = std::make_shared< fastahack::FastaReference >();
		m_fasta_reference->open(this->m_fasta_path);
		setSequence(region);
	}

	FastaReference::~FastaReference()
	{
	}

	void FastaReference::setSequence(Region::SharedPtr region)
	{
		this->m_region = region;

		std::string seqName = this->m_region->getReferenceID();
		if (this->m_region->getStartPosition() == 0 || this->m_region->getEndPosition() == 0)
		{
			this->m_sequence = this->m_fasta_reference->getSequence(seqName);
			this->m_region->setStartPosition(1);
			this->m_region->setEndPosition(this->m_region->getStartPosition() + this->m_sequence.size());
		}
		else
		{
			position startPosition = this->m_region->getStartPosition();
			size_t length = this->m_region->getEndPosition() - startPosition;
			this->m_sequence = this->m_fasta_reference->getSubSequence(seqName, startPosition, length);
		}
	}

} // end namespace gwiz
