#ifndef GWIZ_INODE_H
#define GWIZ_INODE_H

#include <memory>

namespace gwiz
{

	class INode
	{
	public:
		typedef std::shared_ptr<INode> SharedPtr;

		virtual ~INode() {}
	private:

		virtual void ReplaceMe() = 0;

	};

} // end namespace gwiz

#endif // GWIZ_INODE_H
