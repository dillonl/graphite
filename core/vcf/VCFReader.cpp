#include "VCFReader.h"
#include "core/util/Utility.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>


namespace graphite
{
	VCFReader::VCFReader(const std::string& filename, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs, Region::SharedPtr regionPtr, VCFWriter::SharedPtr vcfWriter) :
		m_filename(filename),
		m_region_ptr(nullptr),
		m_file_stream_ptr(nullptr),
		m_preloaded_variant(nullptr),
		m_vcf_writer(vcfWriter)
	{
		openFile(); // open the vcf
		processHeader(bamSamplePtrs); // read the header
		setRegion(regionPtr);
	}

	VCFReader::~VCFReader()
	{
		if (this->m_filename.substr(this->m_filename.find_last_of(".") + 1) == "gz")
		{
			std::static_pointer_cast< igzstream >(m_file_stream_ptr)->close();
		}
		else
		{
			std::static_pointer_cast< std::ifstream >(m_file_stream_ptr)->close();
		}
	}

	void VCFReader::openFile()
	{
		if (this->m_filename.substr(this->m_filename.find_last_of(".") + 1) == "gz")
		{
			auto igzstreamPtr = std::make_shared< igzstream >();
			igzstreamPtr->open(this->m_filename.c_str());
			this->m_file_stream_ptr = igzstreamPtr;
		}
		else
		{
			this->m_file_stream_ptr = std::make_shared< std::ifstream >(this->m_filename, std::fstream::in);
		}
	}

	void VCFReader::setRegion(Region::SharedPtr regionPtr)
	{
		if (regionPtr == nullptr)
		{
			return;
		}
		this->m_region_ptr = regionPtr;
		std::string nextLine;
		while (this->m_preloaded_variant != nullptr)
		{
			// as soon as we are inside the region then break out
			if (this->m_region_ptr->getReferenceID() == this->m_preloaded_variant->getChromosome() &&
				  this->m_region_ptr->getStartPosition() <= this->m_preloaded_variant->getPosition() &&
				  this->m_preloaded_variant->getPosition() <= this->m_region_ptr->getEndPosition())
			{
				break;
			}
			if (!getNextLine(nextLine)) // if we are at the eof then set the variant to nullptr
			{
				this->m_preloaded_variant = nullptr;
				break;
			}
			this->m_preloaded_variant = std::make_shared< Variant >(nextLine, this->m_vcf_writer);
		}
	}

	bool VCFReader::getNextVariants(std::vector< Variant::SharedPtr >& variantPtrs, uint32_t spacing)
	{
		variantPtrs.clear();
		std::string nextLine;
		while (this->m_preloaded_variant != nullptr)
		{
			bool inTheZone = !(this->m_region_ptr != nullptr &&
				!(this->m_region_ptr->getReferenceID() == this->m_preloaded_variant->getChromosome() &&
				  this->m_region_ptr->getStartPosition() <= this->m_preloaded_variant->getPosition() &&
				  this->m_preloaded_variant->getPosition() <= this->m_region_ptr->getEndPosition()));

			bool isNextVariantCloseEnough = (variantPtrs.size() == 0) ? true : ((this->m_preloaded_variant->getChromosome() ==  variantPtrs[variantPtrs.size()-1]->getChromosome()) && (this->m_preloaded_variant->getPosition() -  variantPtrs[variantPtrs.size()-1]->getPosition()) < spacing); // if we are outside the region then break out

			if (!inTheZone || !isNextVariantCloseEnough)
			{
				break;
			}

			variantPtrs.emplace_back(this->m_preloaded_variant);
			if (!getNextLine(nextLine)) // if we are at the eof then set the variant to nullptr
			{
				this->m_preloaded_variant = nullptr;
				break;
			}
			this->m_preloaded_variant = std::make_shared< Variant >(nextLine, this->m_vcf_writer);
		}
		return variantPtrs.size() > 0;
	}

	void VCFReader::processHeader(std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs)
	{
		// for each line read and write it to the VCFWriter
		std::vector< std::string > headerLines;
		std::string line;
		std::string headerColumns = "";
		while (getNextLine(line) && line.c_str()[0] == '#')
		{
			if (line.find("#CHROM") == std::string::npos)
			{
				headerLines.emplace_back(line);
			}
			else
			{
				this->m_vcf_writer->setOriginalVCFSampleNames(line);
				headerColumns = line;
			}
		}
		std::string columnLine = setSamplePtrs(headerColumns, bamSamplePtrs);
		headerLines.emplace_back(columnLine);
		// this->m_vcf_writer->setSamples(columnLine, this->m_sample_ptrs_map);
		this->m_vcf_writer->writeHeader(headerLines);
		if (line.size() > 0)
		{
			this->m_preloaded_variant = std::make_shared< Variant >(line, m_vcf_writer);
		}
	}

	std::string VCFReader::setSamplePtrs(const std::string& columnLine, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs)
	{
		std::string newLine = "";
		std::vector< std::string > sampleNames;
		std::vector< std::string > columns;
		split(columnLine, '\t', columns);
		std::string lastColumnName = "";
		for (auto column : columns)
		{
			std::string upperColumn;
			std::transform(column.begin(), column.end(), std::back_inserter(upperColumn), ::toupper);
			if (std::find(STANDARD_VCF_COLUMN_NAMES.begin(), STANDARD_VCF_COLUMN_NAMES.end(), upperColumn) == STANDARD_VCF_COLUMN_NAMES.end())
			{
				Sample::SharedPtr samplePtr = std::make_shared< Sample >(column, "", "");
				for (auto bamSamplePtr : bamSamplePtrs)
				{
					if (column == bamSamplePtr->getName())
					{
						samplePtr = bamSamplePtr;
						break;
					}
				}
				this->m_sample_ptrs_map.emplace(samplePtr->getName(), samplePtr);
			}
			newLine += column + "\t";
			if (column.size() > 0)
			{
				lastColumnName = column;
			}
		}
		if (lastColumnName.compare("INFO") == 0)
		{
			newLine += "FORMAT\t";
		}

		// add in the rest of the bam samples that aren't in the header already
		for (auto bamSamplePtr : bamSamplePtrs)
		{
			if (this->m_sample_ptrs_map.find(bamSamplePtr->getName()) == this->m_sample_ptrs_map.end())
			{
				this->m_sample_ptrs_map.emplace(bamSamplePtr->getName(), bamSamplePtr);
				newLine += bamSamplePtr->getName() + "\t";
			}
		}
		newLine.pop_back(); // remove the final \t
		return newLine;
	}

}
