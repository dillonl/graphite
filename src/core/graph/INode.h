#ifndef GWIZ_INODE_H
#define GWIZ_INODE_H

#include <memory>

namespace gwiz
{

	class INode
	{
	public:
		typedef std::shared_ptr<INode> SharedPtr;
		INode() {}
	    INode(char* sequence, uint32_t length) :
		    m_sequence(sequence),
			m_length(length)
			{}
		virtual ~INode() {}

	protected:
		char* m_sequence;
		uint32_t m_length;

	};

} // end namespace gwiz

#endif // GWIZ_INODE_H
