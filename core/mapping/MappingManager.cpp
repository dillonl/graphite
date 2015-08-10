#include "MappingManager.h"
#include "IMapping.h"

#include "core/util/ThreadPool.hpp"

#include <iostream>

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
		std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
		auto alignmentPtr = mappingPtr->getAlignmentPtr();
		// std::cout << alignmentPtr->getPosition() << std::endl;
		auto iter = this->m_alignment_mapping_map.find(alignmentPtr);
		if (iter == this->m_alignment_mapping_map.end() || iter->second->getMappingScore() < mappingPtr->getMappingScore())
		{
			this->m_alignment_mapping_map.emplace(alignmentPtr, mappingPtr);
		}
	}

	void MappingManager::evaluateAlignmentMappings(IAdjudicator::SharedPtr adjudicatorPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
		for (auto& iter : this->m_alignment_mapping_map)
		{
			auto funct = std::bind(&IAdjudicator::adjudicateMapping, adjudicatorPtr, iter.second);
			ThreadPool::Instance()->enqueue(funct);
		}
		gwiz::ThreadPool::Instance()->joinAll();
	}
}
