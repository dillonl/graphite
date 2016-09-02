#ifndef GRAPHITE_SEQUENCE_MANAGER_H
#define GRAPHITE_SEQUENCE_MANAGER_H

#include "Sequence.h"

#include <unordered_map>
#include <cstring>
#include <mutex>

namespace graphite
{
	class SequenceManager : private Noncopyable
	{
	public:
		static SequenceManager* Instance();

		inline Sequence::SharedPtr getSequence(const std::string& sequence)
		{
			std::lock_guard< std::mutex > lock(this->m_seq_map_mutex);
			auto sequencePtrIter = this->m_sequence_ptr_map.find(sequence);
			if (sequencePtrIter == this->m_sequence_ptr_map.end())
			{
				auto sequencePtr = std::make_shared< Sequence >(sequence);
				this->m_sequence_ptr_map[sequence] = sequencePtr;
				return sequencePtr;
			}
			else
			{
				return sequencePtrIter->second;
			}
		}

		// created for testing
		size_t getSequenceCount()
		{
			std::lock_guard< std::mutex > lock(this->m_seq_map_mutex);
			return this->m_sequence_ptr_map.size();
		}
		void clearSequences()
		{
			std::lock_guard< std::mutex > lock(this->m_seq_map_mutex);
			this->m_sequence_ptr_map.clear();
		}

		void printSize()
		{
			std::lock_guard< std::mutex > lock(this->m_seq_map_mutex);
			std::cout << "SM size: " << m_sequence_ptr_map.size() << std::endl;
		}

	private:


		SequenceManager() {}
		~SequenceManager() {}

		// contains all the sequences with the string representation as a key
		/* std::unordered_map< const char*, Sequence::SharedPtr > m_sequence_ptr_map; */
		std::mutex m_seq_map_mutex;
		std::unordered_map< std::string, Sequence::SharedPtr > m_sequence_ptr_map;
		std::vector< Sequence::SharedPtr > m_seq_vec;
		static SequenceManager* s_instance;
	};
}

#endif //GRAPHITE_SEQUENCE_MANAGER_H
