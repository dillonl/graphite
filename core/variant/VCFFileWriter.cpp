#include "VCFFileWriter.h"
#include "VCFFileReader.h"
#include "core/file/BGZFFileWriter.h"
#include "core/file/ASCIIFileWriter.h"

namespace graphite
{
	VCFFileWriter::VCFFileWriter(const std::string& filePath, const std::string& originalVCFPath, FileType fileType) :
		m_file_path(filePath), m_header_written(false)
	{
		if (fileType == graphite::FileType::BGZF)
		{
			m_file_writer = std::make_shared< BGZFFileWriter >(filePath);
		}
		else
		{
			m_file_writer = std::make_shared< ASCIIFileWriter >(filePath);
		}
		m_file_writer->open();
	}

	VCFFileWriter::~VCFFileWriter()
	{
	}

	void VCFFileWriter::close()
	{
		m_file_writer->close();
	}

	void VCFFileWriter::writeVariantList(VariantList::SharedPtr variantListPtr, VCFHeader::SharedPtr vcfHeaderPtr)
	{
		if (!m_header_written)
		{
			auto headerStr = vcfHeaderPtr->getHeader();
			m_file_writer->write(headerStr.c_str(), headerStr.size());
			m_header_written = true;
		}
		for(const auto variantPtr : variantListPtr->getAllVariantPtrs())
		{
			auto line = variantPtr->getVariantLine(vcfHeaderPtr);
			m_file_writer->write(line.c_str(), line.size());
		}
	}
}
