#include <stdexcept>
#include "boost/xpressive/xpressive.hpp"

#include "Region.h"

namespace gwiz
{
	using namespace boost::xpressive;

	Region::Region(const std::string& region_string) :
		m_region_string(region_string)
	{
		this->m_end_position = 1;
		if (!boost::spirit::qi::parse(region_string.c_str(), region_string.c_str() + region_string.size(), m_region_parser, this->m_reference_id, this->m_start_position, this->m_end_position))
		{
			throw std::invalid_argument("Region format is invalid");
		}
	}

	Region::~Region()
	{
	}
}
