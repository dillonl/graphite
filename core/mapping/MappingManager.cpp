#include "MappingManager.h"
#include "IMapping.h"

#include "core/util/ThreadPool.hpp"

namespace gwiz
{
	MappingManager* MappingManager::s_instance = nullptr;

	MappingManager::MappingManager()
	{
	}

	MappingManager::~MappingManager()
	{
	}

	MappingManager* MappingManager::Instance()
	{
		if (s_instance == nullptr)
		{
			s_instance = new MappingManager();
		}
		return s_instance;
	}

	void MappingManager::registerMapping(IMapping::SharedPtr mappingPtr)
	{
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		auto iter = this->m_alignment_mapping_map.find(alignmentPtr);
		if (iter == this->m_alignment_mapping_map.end() || iter->second->getMappingScore() < mappingPtr->getMappingScore())
		{
			this->m_alignment_mapping_map.emplace(alignmentPtr, mappingPtr);
		}
	}

	void MappingManager::evaluateAlignmentMappings()
	{
		for (auto& iter : this->m_alignment_mapping_map)
		{
			auto funct = std::bind(&MappingManager::evaluateMapping, this, iter.second);
			ThreadPool::Instance()->enqueue(funct);
		}
	}

	void MappingManager::evaluateMapping(IMapping::SharedPtr mappingPtr)
	{

	}
}
