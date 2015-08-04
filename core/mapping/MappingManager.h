#ifndef GWIZ_MAPPINGMANAGER_H
#define GWIZ_MAPPINGMANAGER_H

#include "IMapping.h"

#include <unordered_map>

namespace gwiz
{
	class MappingManager
	{
	public:
		static MappingManager* Instance();
		MappingManager();
		~MappingManager();

		void registerMapping(IMapping::SharedPtr mappingPtr);
	private:
		MappingManager* s_instance;
		std::unordered_map< IAlignment::SharedPtr, IMapping::SharedPtr > m_alignment_mapping_map;
	};
}

#endif //GWIZ_MAPPINGMANAGER_H
