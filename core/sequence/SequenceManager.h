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

		inline Sequence::SharedPtr getSequence(const char* sequence)
		{
			return this->m_sequence_ptr_map.emplace(sequence, std::make_shared< Sequence >(sequence, strlen(sequence))).first->second;
		}

		// created for testing
		size_t getSequenceCount() { return this->m_sequence_ptr_map.size(); }
		void clearSequences() { this->m_sequence_ptr_map.clear(); }

	private:


		SequenceManager() {}
		~SequenceManager() {}

		// contains all the sequences with the string representation as a key
		std::unordered_map< const char*, Sequence::SharedPtr > m_sequence_ptr_map;
		static SequenceManager* s_instance;
	};
}

#endif //GWIZ_SEQUENCE_MANAGER_H
