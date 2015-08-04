#include "MappingManager.h"
#include "IMapping.h"

namespace gwiz
{
	static MappingManager* s_instance = nullptr;

	MappingManager::MappingManager()
	{
	}

	MappingManager::~MappingManager()
	{
	}

	MappingManager* MappingManager::Instance()
	{
	}

	void MappingManager::registerMapping(IMapping::SharedPtr mappingPtr)
	{
	}
}
