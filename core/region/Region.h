#ifndef GRAPHITE_REGION_H
#define GRAPHITE_REGION_H

#include <boost/noncopyable.hpp>

#include <string>
#include <memory>

#include "core/util/Types.h"
#include "core/parser/RegionParser.hpp"

namespace graphite
{

	/*
	 * This class represents a reference id, start position and an end position.
	 */
	class Region : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< Region > SharedPtr;
		Region(const std::string& regionString);
		Region(const std::string& referenceID, position startPosition, position endPosition);
		~Region();

		std::string getRegionString() const { return m_region_string; }
		std::string getReferenceID() const { return this->m_reference_id; }
		/*
		std::string getReferenceIDNormalized() const
		{
			return this->m_reference_id_normalized;
		}
		*/
		position getStartPosition() const { return this->m_start_position; }
		position getEndPosition() const { return this->m_end_position; }
		void setStartPosition(position startPosition) { this->m_start_position = startPosition; }
		void setEndPosition(position endPosition) { this->m_end_position = endPosition; }

	private:
		std::string m_region_string;
		std::string m_reference_id;
		/* std::string m_reference_id_normalized; */
		position m_start_position;
		position m_end_position;

		RegionParser< const char* > m_region_parser;
	};

}

#endif //GRAPHITE_REGION_H
