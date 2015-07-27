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
		}

		~Sequence(){}

		const char* getSequence() const { return m_seq.c_str(); }
		size_t getLength() const { return m_seq.size(); }

		std::string getSequenceString() const { return m_seq; }
	private:
		std::string m_seq;
	};
}

#endif //GWIZ_SEQUENCE_H
