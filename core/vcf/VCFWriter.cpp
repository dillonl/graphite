#include "VCFWriter.h"
#include "core/util/Utility.h"
#include "core/util/Types.h"

#include <algorithm>

namespace graphite
{
	VCFWriter::VCFWriter(const std::string& filename, const std::string& outputDirectory)
	{
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
		// we need to add the graphite format rows in the header.
		// The code below ensures that we only write the format columns once and in the write place.
		std::vector< std::string > lines;
		bool formatWritten = false;
		for (auto line : headerLines)
		{
			if ((line.find("##FORMAT") != std::string::npos && !formatWritten) || (line.find("#CHROM") != std::string::npos && !formatWritten))
			{
				for (auto formatTuple : this->m_format)
				{
					lines.emplace_back(std::get< 1 >(formatTuple));
				}
				formatWritten = true;
			}
			bool addLine = true;
			for (auto formatTuple : this->m_format) // remove all rows that are in the m_format vector
			{
				if (line.find(std::get< 0 >(formatTuple)) == std::string::npos)
				{
					addLine = false;
					break;
				}
			}
			if (addLine)
			{
				lines.emplace_back(line);
			}
		}
		for (auto line : lines)
		{
			writeLine(line);
		}
		std::string headerLine = "";
		for (auto column : STANDARD_VCF_COLUMN_NAMES)
		{
			headerLine += column + "\t";
		}
		for (uint32_t i = 0; i < this->m_sample_ptrs.size(); ++i)
		{
			headerLine += (i == 0) ? "" : "\t";
			headerLine += this->m_sample_ptrs[i]->getName();
		}
		writeLine(headerLine);
	}

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

	Sample::SharedPtr VCFWriter::getSamplePtr(const std::string& sampleName)
	{
		return this->m_sample_ptrs_map.find(sampleName)->second;
	}

	std::vector< Sample::SharedPtr > VCFWriter::getSamplePtrs()
	{
		return m_sample_ptrs;
	}
}
