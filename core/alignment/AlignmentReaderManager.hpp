#ifndef GRAPHITE_ALIGNMENTREADERMANAGER_HPP
#define GRAPHITE_ALIGNMENTREADERMANAGER_HPP

#include "IAlignmentReader.h"
#include "core/util/Noncopyable.hpp"

#include <string.h>
#include <vector>
#include <memory>
#include <chrono>
#include <thread>
#include <deque>

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
				std::deque< std::shared_ptr< AlignmentReaderType > > alignmentReaderQueue;
				for (auto i = 0; i < numReadersPerFile; ++i)
				{
					auto alignmentReaderPtr = std::make_shared< AlignmentReaderType >(fileName, this);
					alignmentReaderPtr->open();
					alignmentReaderQueue.push_front(alignmentReaderPtr);
				}
				m_alignment_readers_queue_map[fileName] = alignmentReaderQueue;
			}
		}

		~AlignmentReaderManager()
		{
			// close all alignment readers
			for (auto iter : m_alignment_readers_queue_map)
			{
				while (!iter.second.empty())
				{
					iter.second.front()->close();
					iter.second.pop_front();
				}
			}
		}

        typename AlignmentReaderType::SharedPtr getReader(const std::string& fileName)
		{
			auto iter = m_alignment_readers_queue_map.find(fileName);
			if (iter == m_alignment_readers_queue_map.end())
			{
				return nullptr;
			}
			while (true)
			{
				{
					std::lock_guard< std::mutex > readerQueueLockGuard(m_alignment_readers_queue_mutex);
					if (!iter->second.empty())
					{
						auto readerPtr = iter->second.front();
						iter->second.pop_front();
						return readerPtr;
					}
				}
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			return nullptr;
		}


		void checkinReader(typename AlignmentReaderType::SharedPtr readerPtr)
		{
			auto iter = m_alignment_readers_queue_map.find(readerPtr->getPath());
			std::lock_guard< std::mutex > readerQueueLockGuard(m_alignment_readers_queue_mutex);
			iter->second.push_back(readerPtr);
		}

/*

  typedef std::shared_ptr< AlignmentReaderManager > SharedPtr;
  AlignmentReaderManager(std::vector< std::string > fileNames, uint32_t numReadersPerFile) :
  m_current_reader_idx(0)
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
				std::cout << "AlignmentReaderManager says :: make a vector of bools (ischeckedin) and lock the vector and get a non-checked out reader then return it" << std::endl;
				for (auto iter2 : alignmentReaderContainersMap)
				{
					if (iter2.second->m_lock.try_lock())
					{
						std::cout << "reader: " << iter2.second->m_alignment_reader.get() << std::endl;
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
*/

	private:

		// std::unordered_map< std::string, std::unordered_map< uint32_t, std::shared_ptr< AlignmentReaderContainer > > > m_alignment_reader_containers_map;
		std::unordered_map< std::string, std::deque< std::shared_ptr< AlignmentReaderType > > > m_alignment_readers_queue_map;
		std::mutex m_alignment_readers_queue_mutex;
	};
}

#endif //GRAPHITE_ALIGNMENTREADERMANAGER_HPP
