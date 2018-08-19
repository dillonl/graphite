#include "VCFWriter.h"
#include "core/util/Utility.h"
#include "core/util/Types.h"

#include <algorithm>

namespace graphite
{
	VCFWriter::VCFWriter(const std::string& filename, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs, const std::string& outputDirectory) :
		m_bam_sample_ptrs(bamSamplePtrs),
		m_black_format_string(nullptr)
	{
		for (auto samplePtr : m_bam_sample_ptrs)
		{
			m_bam_sample_ptrs_map.emplace(samplePtr->getName(), samplePtr);
		}
		std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
		std::string path(outputDirectory + "/" + base_filename);
		this->m_out_file.open(path);
	}

	VCFWriter::~VCFWriter()
	{
		this->m_out_file.close();
	}

	void VCFWriter::writeLine(const std::string& line)
	{
		this->m_out_file << line.c_str() << std::endl;
	}

	void VCFWriter::writeHeader(const std::vector< std::string >& headerLines)
	{
		m_vcf_column_names.clear();
		m_sample_names.clear();
		// we need to add the graphite format rows in the header.
		// The code below ensures that we only write the format columns once and in the write place.
		std::vector< std::string > lines;
		bool formatWritten = false;
		for (auto line : headerLines)
		{
			bool isHeaderCols = false;
			if (line.find("#CHROM") != std::string::npos)
			{
				split(line, '\t', this->m_vcf_column_names);
				isHeaderCols = true;
			}
			if ((line.find("##FORMAT") != std::string::npos && !formatWritten) || (line.find("#CHROM") != std::string::npos && !formatWritten))
			{
				for (auto formatTuple : this->m_format)
				{
					lines.emplace_back(std::get< 1 >(formatTuple));
				}
				formatWritten = true;
			}
			if (!isHeaderCols)
			{
				lines.emplace_back(line);
			}
		}
		std::unordered_set< std::string > writtenLines;
		for (auto line : lines)
		{
			if (writtenLines.find(line) == writtenLines.end())
			{
				writeLine(line);
				writtenLines.emplace(line);
			}
		}
		std::string headerLine = "";
		bool first = true;
		for (auto headerName : this->m_vcf_column_names)
		{
			if (!first)
			{
				headerLine += "\t";
			}
			headerLine += headerName;
			first = false;
			if (STANDARD_VCF_COLUMN_NAMES_SET.find(headerName) == STANDARD_VCF_COLUMN_NAMES_SET.end())
			{
				this->m_sample_names.emplace(this->m_sample_names.end(), headerName);
			}
		}
		first = true;
		for (auto bamSamplePtr : this->m_bam_sample_ptrs)
		{
			auto sampleName = bamSamplePtr->getName();
			auto iter = std::find_if(this->m_vcf_column_names.begin(), this->m_vcf_column_names.end(), [&sampleName](const std::string& columnName)
									 {
										 return columnName.compare(sampleName) == 0;
									 });
			if (iter == this->m_vcf_column_names.end()) // if sample not in vcf
			{
				if (!first)
				{
					headerLine += "\t";
				}
				headerLine += sampleName;
				this->m_vcf_column_names.emplace(this->m_vcf_column_names.end(), sampleName);
				this->m_sample_names.emplace(this->m_sample_names.end(), sampleName);
			}
			else
			{
				this->m_sample_name_in_vcf.emplace(sampleName, true);
			}
			first = false;
		}
		writeLine(headerLine);
	}

	/*
	void VCFWriter::setSamples(const std::string& columnHeaderLine, std::unordered_map< std::string, Sample::SharedPtr >& samplePtrsMap)
	{
		this->m_sample_ptrs.clear();
		this->m_sample_ptrs_map.clear();
		std::vector< std::string > columns;
		split(columnHeaderLine, '\t', columns);
		for (auto column : columns)
		{
			std::string upperColumn;
			std::transform(column.begin(), column.end(), std::back_inserter(upperColumn), ::toupper);
			if (std::find(STANDARD_VCF_COLUMN_NAMES.begin(), STANDARD_VCF_COLUMN_NAMES.end(), upperColumn) == STANDARD_VCF_COLUMN_NAMES.end())
			{
				auto iter = samplePtrsMap.find(column);
				if (iter == samplePtrsMap.end())
				{
					std::cout << "error in VCFWriter::setSamples sample: " << column << " not found" << std::endl;
				}
				else
				{
					this->m_sample_ptrs.emplace_back(iter->second);
					this->m_sample_ptrs_map.emplace(iter->first, iter->second);
				}
			}
		}
	}
	*/

	Sample::SharedPtr VCFWriter::getSamplePtr(const std::string& sampleName)
	{
		return this->m_bam_sample_ptrs_map.find(sampleName)->second;
	}

	std::vector< std::string > VCFWriter::getSampleNames()
	{
		return m_sample_names;
	}

	bool VCFWriter::isSampleNameInOriginalVCF(const std::string& sampleName)
	{
		return (m_sample_name_in_vcf.find(sampleName) != m_sample_name_in_vcf.end());
	}

	bool VCFWriter::isSampleNameInBam(const std::string& sampleName)
	{
		return (m_bam_sample_ptrs_map.find(sampleName) != m_bam_sample_ptrs_map.end());
	}

	std::vector< std::string > VCFWriter::getColumnNames()
	{
		return this->m_vcf_column_names;
	}

	void VCFWriter::setBlankFormatString(const std::string& blankFormatString)
	{
		if (this->m_black_format_string == nullptr)
		{
			this->m_black_format_string = std::make_shared< std::string >(blankFormatString);
		}
	}

	std::shared_ptr< std::string > VCFWriter::getBlankFormatStringPtr()
	{
		return this->m_black_format_string;
	}
}
