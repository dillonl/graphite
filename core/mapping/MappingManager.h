#ifndef GWIZ_MAPPINGMANAGER_H
#define GWIZ_MAPPINGMANAGER_H

#include "IMapping.h"

#include <unordered_map>

#include <boost/noncopyable.hpp>

namespace gwiz
{
	class MappingManager : private boost::noncopyable
	{
	public:
		static MappingManager* Instance();

		/*
		 * Checks the mapping's mappingScore and only
		 * adds the alignmentPtr if the passed in
		 * mappingPtr's mapping score is larger.
		 */
		void registerMapping(IMapping::SharedPtr mappingPtr);
		void evaluateAlignmentMappings();
	private:
		void evaluateMapping(IMapping::SharedPtr mappingPtr);
		MappingManager();
		~MappingManager();

		static MappingManager* s_instance;
		std::unordered_map< IAlignment::SharedPtr, IMapping::SharedPtr > m_alignment_mapping_map;
	};
}

#endif //GWIZ_MAPPINGMANAGER_H
