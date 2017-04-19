#include "FastaReference.h"


namespace graphite
{

	FastaReference::FastaReference(const std::string& path, Region::SharedPtr region) :
		m_fasta_path(path),
		IReference()
	{
		m_fasta_reference = std::make_shared< ::FastaReference >();
		m_fasta_reference->open(this->m_fasta_path);
		setSequence(region);
	}

	FastaReference::~FastaReference()
	{
	}

	void FastaReference::setSequence(Region::SharedPtr regionPtr)
	{
		std::string seqName = regionPtr->getReferenceID();
		this->m_sequence = this->m_fasta_reference->getSequence(seqName);

		this->m_region = std::make_shared< Region >(seqName, Region::BASED::ONE);
		this->m_region->setBased(Region::BASED::ZERO);
		this->m_region->setStartPosition(0);
		this->m_region->setEndPosition(this->m_region->getStartPosition() + this->m_sequence.size());

		/*
		if (this->m_region->getStartPosition() == 0 || this->m_region->getEndPosition() == 0)
		{
			this->m_sequence = this->m_fasta_reference->getSequence(seqName);
			this->m_region->setStartPosition(0);
			this->m_region->setEndPosition(this->m_region->getStartPosition() + this->m_sequence.size());
		}
		else
		{
			if (regionPtr->getBased() == Region::BASED::ONE) // fastahack is one based so we just provide the start position as is, length we must add one because of math!
			{
				this->m_sequence = this->m_fasta_reference->getSubSequence(seqName, this->m_region->getStartPosition(), (this->m_region->getEndPosition() - this->m_region->getStartPosition()) + 1);
				this->m_region->setStartPosition(this->m_region->getStartPosition() - 1);
				this->m_region->setEndPosition(this->m_region->getEndPosition() - 1);
			}
			else
			{
				this->m_sequence = this->m_fasta_reference->getSubSequence(seqName, this->m_region->getStartPosition() + 1, (this->m_region->getEndPosition() - this->m_region->getStartPosition()) + 1);
			}
		}
		*/
	}

} // end namespace graphite
