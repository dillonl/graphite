#include <stdexcept>
#include "boost/xpressive/xpressive.hpp"

#include "Region.h"

namespace gwiz
{
	using namespace boost::xpressive;

	Region::Region(const std::string& region_string) :
		m_region_string(region_string)
	{
		sregex rex = sregex::compile("^(\\w+):(\\d+)-(\\d+)$"); // in the form of 20:1000-2000

		smatch matches;

		if(regex_match(region_string, matches, rex))
		{
			this->m_reference_id = matches[1];
			this->m_start_position = stoi(matches[2]);
			this->m_end_position = stoi(matches[3]);
		}
		else
		{
			throw std::invalid_argument("Region string is invalid");
		}
	}

	Region::~Region()
	{
	}

}
