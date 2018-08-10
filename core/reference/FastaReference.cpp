#include "FastaReference.h"

namespace graphite
{
	FastaReference::FastaReference(const std::string& path) :
		m_current_reference_id(""),
		m_reference_sequence("")
	{
		m_fasta_reference = std::make_shared< ::FastaReference >();
		m_fasta_reference->open(path);
	}

	FastaReference::~FastaReference()
	{
		// m_fasta_reference is closed when its destructor is called so we don't worry about it
	}

	void FastaReference::setReferenceIDAndSequence(Region::SharedPtr regionPtr)
	{
		std::string referenceID = regionPtr->getReferenceID();
		if (strcmp(referenceID.c_str(), this->m_current_reference_id.c_str()) != 0)
		{
			this->m_current_reference_id = referenceID;
			this->m_reference_sequence = this->m_fasta_reference->getSequence(referenceID);
		}
	}

	std::string FastaReference::getSequenceStringFromRegion(Region::SharedPtr regionPtr)
	{
		setReferenceIDAndSequence(regionPtr);
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();
		if (regionPtr->getBased() == Region::BASED::ONE)
		{
			startPosition -= 1;
			endPosition -= 1;
		}
		position length = endPosition - startPosition;
		return std::string(this->m_reference_sequence.c_str() + startPosition, length);
	}

	const char* FastaReference::getSequenceFromRegion(Region::SharedPtr regionPtr)
	{
		setReferenceIDAndSequence(regionPtr);
		position startPosition = regionPtr->getStartPosition();
		if (regionPtr->getBased() == Region::BASED::ONE)
		{
			startPosition -= 1;
		}
		return &this->m_reference_sequence.c_str()[startPosition];
	}
}
