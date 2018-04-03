#ifndef GRAPHITE_VCFFILEWRITER_H
#define GRAPHITE_VCFFILEWRITER_H

#include "core/util/Noncopyable.hpp"
#include "core/file/IFileWriter.h"
#include "core/file/IFile.h"
#include "core/variant/VariantList.h"

#include <string>

namespace graphite
{
    class VCFFileWriter : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VCFFileWriter > SharedPtr;
        VCFFileWriter(const std::string& filePath, const std::string& originalVCFPath, FileType fileType);
		~VCFFileWriter();

		void writeVariantList(VariantList::SharedPtr variantListPtr, VCFHeader::SharedPtr vcfHeader);
		void close();

	private:
		IFileWriter::SharedPtr m_file_writer;
		std::string m_file_path;
		bool m_header_written;

	};
}

#endif //GRAPHITE_VCFFILEWRITER_H
