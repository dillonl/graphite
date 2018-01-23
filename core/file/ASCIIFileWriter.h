#ifndef GRAPHITE_ASCIIFILEWRITER
#define GRAPHITE_ASCIIFILEWRITER

#include <iostream>
#include <fstream>

#include "IFileWriter.h"

namespace graphite
{
	class ASCIIFileWriter : public IFileWriter
	{
	public:
		typedef std::shared_ptr< ASCIIFileWriter > SharedPtr;
		ASCIIFileWriter(const std::string& path);
		~ASCIIFileWriter() override;

		bool open() override;
		void close() override;
		bool write(const char* data, size_t dataLength) override;
		std::string getFilePath() override { return m_file_path; }

	protected:
		std::ofstream m_out_stream;
		std::string m_file_path;
		bool m_opened;
	};
}
#endif //GRAPHITE_ASCIIFILEWRITER
