#ifndef GRAPHITE_INODE_H
#define GRAPHITE_INODE_H

#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"

#include <string>
#include <memory>
#include <vector>

namespace graphite
{

	class INode : private Noncopyable
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
		size_t getLength() { return m_length; }
		uint32_t getID() { return m_id; }
		position getPosition() { return m_position; }

	protected:
		const char* m_sequence;
		size_t m_length;
		uint32_t m_id;
		position m_position;
		std::vector< INode::SharedPtr > m_previous_nodes;
		std::vector< INode::SharedPtr > m_next_nodes;
	};

} // end namespace graphite

#endif // GRAPHITE_INODE_H
