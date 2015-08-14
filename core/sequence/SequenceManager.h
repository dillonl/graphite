#ifndef GRAPHITE_SEQUENCE_MANAGER_H
#define GRAPHITE_SEQUENCE_MANAGER_H

#include "Sequence.h"

#include <boost/noncopyable.hpp>
#include <unordered_map>
#include <cstring>
#include <mutex>

namespace graphite
{
	class SequenceManager : private boost::noncopyable
	{
	public:
		static SequenceManager* Instance();

		inline Sequence::SharedPtr getSequence(const std::string& sequence)
		{
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

#endif //GRAPHITE_SEQUENCE_MANAGER_H
