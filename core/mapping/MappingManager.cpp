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
		// adjudicateMappings(adjudicatorPtr);
		incrementVariantCounts();
	}

	void MappingManager::adjudicateMappings(IAdjudicator::SharedPtr adjudicatorPtr)
	{
		std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
		std::deque< std::shared_ptr< std::future< void > > > futureFuncts;
		for (auto& mappingPtr : this->m_mappings)
		{
			// std::cout << "mapping" << std::endl;
			auto funct = std::bind(&IAdjudicator::adjudicateMapping, adjudicatorPtr, mappingPtr);
			futureFuncts.push_back(ThreadPool::Instance()->enqueue(funct));
		}
		while (!futureFuncts.empty())
		{
			futureFuncts.front()->wait();
			futureFuncts.pop_front();
		}
	}

	void MappingManager::incrementVariantCounts()
	{
		std::lock_guard< std::mutex > lock(this->m_alignment_mapping_map_mutex);
		std::deque< std::shared_ptr< std::future< void > > > futureFuncts;
		for (auto& mappingPtr : this->m_mappings)
		{
			auto funct = std::bind(&IMapping::incrementAlleleCounts, mappingPtr);
			futureFuncts.push_back(ThreadPool::Instance()->enqueue(funct));
			/*
			if (mappingPtr->getMapped())
			{
				mappingPtr->printMapping();
			}
			*/
		}
		while (!futureFuncts.empty())
		{
			futureFuncts.front()->wait();
			futureFuncts.pop_front();
		}
	}

	void MappingManager::clearRegisteredMappings()
	{
		this->m_mappings.clear();
		this->m_alignment_mapping_map.clear();
	}
}
