#ifndef GRAPHITE_REGION_H
#define GRAPHITE_REGION_H

#include <string>
#include <memory>

#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"

/* #include "core/parser/RegionParser.hpp" */

namespace graphite
{

	/*
	 * This class represents a reference id, start position and an end position.
	 */
	class Region : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Region > SharedPtr;

		enum class BASED { ZERO = 0, ONE = 1 }; // based on this https://www.biostars.org/p/84686/ interpretation of zero and one based

		Region(const std::string& regionString, BASED based);
		Region(const std::string& referenceID, position startPosition, position endPosition, BASED based);
		~Region();

		std::string getRegionString() const { return m_region_string; }
		std::string getReferenceID() const { return this->m_reference_id; }
		position getStartPosition() const { return this->m_start_position; }
		position getEndPosition() const { return this->m_end_position; }
		void setStartPosition(position startPosition) { this->m_start_position = startPosition; }
		void setEndPosition(position endPosition) { this->m_end_position = endPosition; }
		BASED getBased() { return m_based; }
		void setBased(BASED based);

	private:
		std::string m_region_string;
		std::string m_reference_id;
		position m_start_position;
		position m_end_position;
		BASED m_based;

		/* RegionParser< const char* > m_region_parser; */
	};

}

#endif //GRAPHITE_REGION_H
