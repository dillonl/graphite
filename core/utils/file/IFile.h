#ifndef GWIZ_IFILE_H
#define GWIZ_IFILE_H

#include <boost/noncopyable.hpp>
#include <string>
#include <memory>

namespace gwiz
{
	class IFile : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr<IFile> SharedPtr;
	    IFile(const std::string& path)
			: m_file_path(path), m_opened(false), m_file_size(0)
		{
		}
		virtual ~IFile(){}

		virtual void Open() = 0;
		virtual void Close() = 0;

	protected:
		inline bool fileExists(const std::string& name)
		{
		    if (FILE *file = fopen(name.c_str(), "r"))
			{
				fclose(file);
				return true;
			}
			else
			{
				return false;
			}
		}

		std::string m_file_path;
		bool m_opened;
		size_t m_file_size;
	};
} // end namespace gwiz

#endif // GWIZ_IFILE_H
