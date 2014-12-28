#ifndef GWIZ_REGION_H
#define GWIZ_REGION_H

#include <string>

#include "utils/NonCopyable.h"
#include "utils/Types.h"

namespace gwiz
{

	/*
	 * This class represents a reference id, start position and an end position.
	 */
	class Region : private noncopyable
	{
	public:
		typedef std::shared_ptr<Region> SharedPtr;
		Region(const std::string& region_string);
		Region() = delete; // force region_string usage
		~Region();

		std::string getRegionString() const { return m_region_string; }
		std::string getReferenceID() const { return this->m_reference_id; }
		position getStartPosition() const { return this->m_start_position; }
		position getEndPosition() const { return this->m_end_position; }

	private:
		std::string m_region_string;
		std::string m_reference_id;
		position m_start_position;
		position m_end_position;
	};

}

#endif //GWIZ_REGION_H
