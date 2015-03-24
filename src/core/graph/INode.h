#ifndef GWIZ_INODE_H
#define GWIZ_INODE_H

#include <boost/noncopyable.hpp>

#include "core/utils/Types.h"
#include <string>
#include <memory>

namespace gwiz
{

	class INode : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr<INode> SharedPtr;
		INode() {}
    	INode(const char* sequence, size_t length) :
		        m_sequence(sequence),
				m_length(length)
		{
		}

		virtual ~INode() {}

		virtual const char* getSequence() { return m_sequence; }
		uint32_t getLength() { return m_length; }
		position getPosition() { return m_position; }

	protected:
		const char* m_sequence;
		size_t m_length;
		position m_position;
	};

} // end namespace gwiz

#endif // GWIZ_INODE_H
