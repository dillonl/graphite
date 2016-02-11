#include "MappingManager.h"
#include "IMapping.h"

#include "core/util/ThreadPool.hpp"

#include <iostream>

namespace graphite
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
		this->m_mappings.emplace_back(mappingPtr);
	}

	void MappingManager::evaluateAlignmentMappings(IAdjudicator::SharedPtr adjudicatorPtr)
	{
		{
			std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
			for (auto& mappingPtr : this->m_mappings)
			{
				auto funct = std::bind(&IAdjudicator::adjudicateMapping, adjudicatorPtr, mappingPtr);
				ThreadPool::Instance()->enqueue(funct);
			}
		}
		ThreadPool::Instance()->joinAll();
		{
			std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
			for (auto& mappingPtr : this->m_mappings)
			{
				mappingPtr->incrementAlleleCounts();
				if (mappingPtr->getMapped())
				{
					mappingPtr->printMapping();
				}
			}
		}
	}

	void MappingManager::clearRegisteredMappings()
	{
		this->m_mappings.clear();
		this->m_alignment_mapping_map.clear();
	}
}
