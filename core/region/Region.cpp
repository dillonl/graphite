#include "Region.h"

#include "core/util/Utility.h"

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

namespace graphite
{

	Region::Region(const std::string& regionString) :
		m_region_string(regionString),
		m_start_position(0),
		m_end_position(0)
	{
		std::vector< std::string > chromWithPositionComponents;
		std::vector< std::string > positionComponents;
		split(regionString, ':', chromWithPositionComponents);

		if (chromWithPositionComponents.size() == 0)
		{
			this->m_reference_id = regionString;
		}
		else
		{
			this->m_reference_id = chromWithPositionComponents[0];
			split(chromWithPositionComponents[1], '-', positionComponents);

			if (positionComponents.size() == 2)
			{
				m_start_position = std::stol(positionComponents[0]);
				m_end_position = std::stol(positionComponents[1]);
			}
		}
		if (this->m_start_position > this->m_end_position || regionString.size() == 0)
		{
			throw std::invalid_argument("Region format is invalid");
		}
		if (this->m_start_position == 0 && this->m_end_position == 0)
		{
			this->m_end_position = MAX_POSITION;
		}
	}


	Region::Region(const std::string& referenceID, position startPosition, position endPosition) :
	    m_reference_id(referenceID),
		m_start_position(startPosition),
		m_end_position(endPosition),
		m_region_string(referenceID + ":" + std::to_string(startPosition) + "-" + std::to_string(endPosition))
	{
	}

	Region::~Region()
	{
	}
}
