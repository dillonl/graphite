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
			// we must subtract 1 because fastahack is 1 based and strings are 0 based and since we
			// are using this as a char* then we must account for that
			position startPosition = this->m_region->getStartPosition() - 1;
			// we must add 1 because fastahack is not inclusive
			// so [100-200] would be position 100 to 199
			size_t length = (this->m_region->getEndPosition() + 1) - startPosition;
			this->m_sequence = this->m_fasta_reference->getSubSequence(seqName, startPosition, length);
		}
	}

} // end namespace gwiz
