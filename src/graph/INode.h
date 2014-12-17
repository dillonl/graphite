#ifndef GWIZ_INODE_H
#define GWIZ_INODE_H

namespace gwiz
{

	class INode
	{
	public:
		virtual ~INode() {}

	private:

		virtual void ReplaceMe() = 0;

	};

} // end namespace gwiz

#endif // GWIZ_INODE_H
