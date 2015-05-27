#ifndef GWIZ_REGION_H
#define GWIZ_REGION_H

#include <boost/noncopyable.hpp>

#include <string>
#include <memory>

#include "core/util/Types.h"
#include "RegionParser.hpp"

namespace gwiz
{

	/*
	 * This class represents a reference id, start position and an end position.
	 */
	class Region : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< Region > SharedPtr;
		Region(const std::string& regionString);
		Region() = delete; // force region_string usage
		~Region();

		std::string getRegionString() const { return m_region_string; }
		std::string getReferenceID() const { return this->m_reference_id; }
		position getStartPosition() const { return this->m_start_position; }
		position getEndPosition() const { return this->m_end_position; }

		void setStartPosition(position startPosition) { this->m_start_position = startPosition; }
		void setEndPosition(position endPosition) { this->m_end_position = endPosition; }

	private:
		std::string m_region_string;
		std::string m_reference_id;
		position m_start_position;
		position m_end_position;

		RegionParser< const char* > m_region_parser;
	};

}

#endif //GWIZ_REGION_H
