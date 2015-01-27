#ifndef GWIZ_INODE_H
#define GWIZ_INODE_H

#include "core/utils/Types.h"
#include <memory>

namespace gwiz
{

	class INode
	{
	public:
		typedef std::shared_ptr<INode> SharedPtr;
		INode() {}
    	INode(const char* sequence, size_t length) :
		    m_sequence(sequence),
			m_length(length)
		{
			nodeSeq = std::string(m_sequence, m_length);
		}
			/*
    	INode(const char* sequence, position pos,  uint32_t length) :
		    m_sequence(sequence),
				m_position(pos),
				m_length(length)
			{
				setSequence();
			}
			*/
		virtual ~INode() {}

		const char* getSequence() { return m_sequence; }
		uint32_t getLength() { return m_length; }

		void setLength(uint32_t length) { m_length = length; setSequence(); }

		std::string nodeSeq;
		std::string sequence;
		const char* test;

	protected:
		void setSequence()
		{
			sequence = std::string(m_sequence, m_length);
		}

		const char* m_sequence;
		size_t m_length;
		position m_position;

	};

} // end namespace gwiz

#endif // GWIZ_INODE_H
