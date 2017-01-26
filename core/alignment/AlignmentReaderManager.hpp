#ifndef GRAPHITE_ALIGNMENTREADERMANAGER_HPP
#define GRAPHITE_ALIGNMENTREADERMANAGER_HPP

#include "IAlignmentReader.h"
#include "core/util/Noncopyable.hpp"

#include <string.h>
#include <vector>
#include <memory>
#include <chrono>
#include <thread>

namespace graphite
{
	template < class AlignmentReaderType >
	class AlignmentReaderManager : private Noncopyable
	{
	private:
		struct AlignmentReaderContainer
		{
			std::mutex m_lock;
			typename AlignmentReaderType::SharedPtr m_alignment_reader;
		};
	public:
		typedef std::shared_ptr< AlignmentReaderManager > SharedPtr;
		AlignmentReaderManager(std::vector< std::string > fileNames, uint32_t numReadersPerFile)
		{
			for (auto& fileName : fileNames)
			{
				std::unordered_map< uint32_t, std::shared_ptr< AlignmentReaderContainer > > alignmentReaderContainers;
				for (auto i = 0; i < numReadersPerFile; ++i)
				{
					auto alignmentContainerPtr = std::make_shared< AlignmentReaderContainer >();
					auto alignmentReaderPtr = std::make_shared< AlignmentReaderType >(fileName, this);
					alignmentReaderPtr->open();
					alignmentContainerPtr->m_alignment_reader = alignmentReaderPtr;
					alignmentReaderContainers[alignmentReaderPtr->getReaderID()] = alignmentContainerPtr;
				}
				m_alignment_reader_containers_map[fileName] = alignmentReaderContainers;
			}
		}

		~AlignmentReaderManager()
		{
		}

		typename AlignmentReaderType::SharedPtr getReader(const std::string& fileName)
		{
			auto iter = m_alignment_reader_containers_map.find(fileName);
			if (iter == m_alignment_reader_containers_map.end())
			{
				return nullptr;
			}
			while (true)
			{
				std::unordered_map< uint32_t, std::shared_ptr< AlignmentReaderContainer > > alignmentReaderContainersMap = iter->second;
				for (auto iter2 : alignmentReaderContainersMap)
				{
					if (iter2.second->m_lock.try_lock())
					{
						return iter2.second->m_alignment_reader;
					}
				}
				// std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			return nullptr;
		}

		void checkinReader(typename AlignmentReaderType::SharedPtr readerPtr)
		{
			auto iter = m_alignment_reader_containers_map.find(readerPtr->getPath());
			if (iter != m_alignment_reader_containers_map.end())
			{
				auto readerIter = iter->second.find(readerPtr->getReaderID());
				if (readerIter != iter->second.end())
				{
					readerIter->second->m_lock.unlock();
				}
			}
		}

	private:

		std::unordered_map< std::string, std::unordered_map< uint32_t, std::shared_ptr< AlignmentReaderContainer > > > m_alignment_reader_containers_map;
	};
}

#endif //GRAPHITE_ALIGNMENTREADERMANAGER_HPP
