#ifndef GWIZ_SEQUENCE_H
#define GWIZ_SEQUENCE_H

#include <boost/noncopyable.hpp>

#include <memory>
#include <string>
#include <vector>
#include <iostream>

namespace gwiz
{
	class Sequence : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< Sequence > SharedPtr;

	    Sequence(const std::string& seq) :
		    m_seq(seq)
		{
			static uint32_t s_id1 = 0;
			id = s_id1++;
		}

	    Sequence(const char* seq, size_t len) :
			m_seq(std::string(seq, len))
		{
			static uint32_t s_id2 = 1000000;
			id = s_id2++;
		}

		~Sequence()
		{
			std::cout << "deleting sequence: " << id << std::endl;
		}

		const char* getSequence() const { return m_seq.c_str(); }
		size_t getLength() const { return m_seq.size(); }

		std::string getSequenceString() const { return m_seq; }

		uint32_t id;
	private:
		std::string m_seq;
	};
}

#endif //GWIZ_SEQUENCE_H
