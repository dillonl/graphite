#ifndef GWIZ_SEQUENCE_H
#define GWIZ_SEQUENCE_H

#include <boost/noncopyable.hpp>

#include <memory>
#include <string>
#include <vector>

namespace gwiz
{
	class Sequence : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< Sequence > SharedPtr;
	    Sequence(const char* seq, size_t seqLen) :
		    m_seq(seq, seqLen)
		{
		}

		~Sequence() {}

		const char* getSequence() const { return m_seq.c_str(); }
		size_t getLength() const { return m_seq.size(); }
	private:
		std::string m_seq;
	};
}

#endif //GWIZ_SEQUENCE_H
