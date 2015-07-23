#ifndef GWIZ_SEQUENCE_MANAGER_H
#define GWIZ_SEQUENCE_MANAGER_H

#include "Sequence.h"

#include <boost/noncopyable.hpp>
#include <unordered_map>
#include <cstring>

namespace gwiz
{
	class SequenceManager : private boost::noncopyable
	{
	public:
		static SequenceManager* Instance();

		inline Sequence::SharedPtr getSequence(const std::string& sequence)
		{
			/* m_seq_vec.push_back(sequencePtr); */
			auto sequencePtrIter = this->m_sequence_ptr_map.find(sequence);
			if (sequencePtrIter == this->m_sequence_ptr_map.end())
			{
				auto sequencePtr = std::make_shared< Sequence >(sequence);
				this->m_sequence_ptr_map[sequence] = sequencePtr;
				std::cout << sequencePtr->id << std::endl;
				return sequencePtr;
			}
			else
			{
				std::cout << sequencePtrIter->second->id << std::endl;
				return sequencePtrIter->second;
			}
		}

		inline Sequence::SharedPtr getSequence(const char* sequence, size_t length)
		{
			/*
			auto sequencePtr = std::make_shared< Sequence >(sequence, length);
			std::cout << sequencePtr->id << std::endl;
			return this->m_sequence_ptr_map.emplace(sequencePtr->getSequenceString(), sequencePtr).first->second;
			*/

			std::cout << "t2" << std::endl;
			std::string sequenceString = std::string(sequence, length);
			auto sequencePtrIter = this->m_sequence_ptr_map.find(sequenceString);
			if (sequencePtrIter == this->m_sequence_ptr_map.end())
			{
				auto sequencePtr = std::make_shared< Sequence >(sequenceString);
				this->m_sequence_ptr_map[sequenceString] = sequencePtr;
				std::cout << sequencePtr->id << " new" << std::endl;
				return sequencePtr;
			}
			else
			{
				std::cout << sequencePtrIter->second->id << " old" << std::endl;
				return sequencePtrIter->second;
			}
		}

		// created for testing
		size_t getSequenceCount() { return this->m_sequence_ptr_map.size(); }
		void clearSequences() { this->m_sequence_ptr_map.clear(); }

		void printSize()
		{
			static std::mutex mtx;
			mtx.lock();
			std::cout << "SM size: " << m_sequence_ptr_map.size() << std::endl;
			mtx.unlock();
		}

	private:


		SequenceManager() {}
		~SequenceManager() {}

		// contains all the sequences with the string representation as a key
		/* std::unordered_map< const char*, Sequence::SharedPtr > m_sequence_ptr_map; */
		std::unordered_map< std::string, Sequence::SharedPtr > m_sequence_ptr_map;
		std::vector< Sequence::SharedPtr > m_seq_vec;
		static SequenceManager* s_instance;
	};
}

#endif //GWIZ_SEQUENCE_MANAGER_H
