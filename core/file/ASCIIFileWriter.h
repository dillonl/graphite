#ifndef GRAPHITE_ASCIIFILEWRITER
#define GRAPHITE_ASCIIFILEWRITER

#include <iostream>
#include <fstream>

#include "IFile.h"

namespace graphite
{
	class ASCIIFileWriter : public IFile
	{
	public:
		typedef std::shared_ptr< ASCIIFileWriter > SharedPtr;
		ASCIIFileWriter(const std::string& path);
		~ASCIIFileWriter() override;

		void Write(const char* output, const size_t size);
		virtual void Open() override;
		void Close() override;

	private:
		std::ofstream m_out_stream;
	};
}
#endif //GRAPHITE_ASCIIFILEWRITER
