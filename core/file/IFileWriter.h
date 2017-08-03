#ifndef GRAPHITE_IFILE_WRITER_H
#define GRAPHITE_IFILE_WRITER_H

#include "core/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	enum class FileType {BGZF, ASCII, SAM};

	class IFileWriter : private Noncopyable
	{
	public:
		typedef std::shared_ptr< IFileWriter > SharedPtr;

	    IFileWriter(FileType ft) : m_file_type(ft) {}
		virtual ~IFileWriter() {}

		virtual bool open() = 0;
		virtual void close() = 0;
		virtual bool write(const char* data, size_t dataLength) = 0;

		FileType getFileType() { return m_file_type; }
		virtual std::string getFilePath() = 0;

	protected:
		FileType m_file_type;

	private:
	};
}

#endif //GRAPHITE_IFILE_WRITER_H
