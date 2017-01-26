#ifndef GRAPHITE_IFILE_H
#define GRAPHITE_IFILE_H

#include <string>
#include <memory>
#include <fstream>

#include "core/util/Noncopyable.hpp"

/* #include <boost/filesystem.hpp> */

namespace graphite
{
	class IFile : private Noncopyable
	{
	public:
		typedef std::shared_ptr<IFile> SharedPtr;
	    IFile(const std::string& path)
			: m_file_path(path), m_opened(false)//, m_file_size(0)
		{
		}
		virtual ~IFile(){}

		virtual void Open() = 0;
		virtual void Close() = 0;
		/* virtual const char* getNextLine() = 0; */
		virtual bool getNextLine(std::string& line) = 0;
		virtual void setFilePosition(uint64_t pos) = 0;

	protected:
		inline bool fileExists(const std::string& filename)
		{
			/* return boost::filesystem::exists(name); */
			std::ifstream ifile(filename.c_str());
			return (bool)ifile;
		}

		std::string m_file_path;
		bool m_opened;
		/* size_t m_file_size; */
	};
} // end namespace graphite

#endif // GRAPHITE_IFILE_H
