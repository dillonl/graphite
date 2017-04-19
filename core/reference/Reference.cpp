#include "Reference.h"

namespace graphite
{

	Reference::Reference(const std::string& regionSequence, Region::SharedPtr regionPtr)
	{
		this->m_region = std::make_shared< Region >(regionPtr->getRegionString(), Region::BASED::ONE);
		this->m_sequence = regionSequence;
		this->m_region->setBased(Region::BASED::ZERO);
		if (regionPtr->getBased() == Region::BASED::ONE)
		{
			this->m_region->setStartPosition(this->m_region->getStartPosition() - 1);
			this->m_region->setEndPosition(this->m_region->getEndPosition() - 1);
		}
	}

	Reference::~Reference()
	{
	}

}
