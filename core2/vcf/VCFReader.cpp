#include "VCFReader.h"
#include "core2/util/Utility.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>


namespace graphite
{
	VCFReader::VCFReader(const std::string& filename, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs, Region::SharedPtr regionPtr) :
		m_filename(filename),
		m_region_ptr(nullptr),
		m_file_stream_ptr(nullptr),
		m_preloaded_variant(nullptr)
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
			// if we are outside the region then break out
			if (this->m_region_ptr != nullptr &&
				!(this->m_region_ptr->getReferenceID() == this->m_preloaded_variant->getChromosome() &&
				  this->m_region_ptr->getStartPosition() <= this->m_preloaded_variant->getPosition() &&
				  this->m_preloaded_variant->getPosition() <= this->m_region_ptr->getEndPosition()))
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

	void VCFReader::registerVCFWriter(VCFWriter::SharedPtr vcfWriter)
	{
		this->m_vcf_writer = vcfWriter;
	}

	void VCFReader::getRegionsFromVCF(std::vector< Region::SharedPtr >& regionPtrs)
	{

	}

	void VCFReader::processHeader(std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs)
	{
		// for each line read and write it to the VCFWriter
		std::vector< std::string > headerLines;
		std::string line;
		std::string headerColumns = "";
		while (getNextLine(line) && line.c_str()[0] == '#')
		{
			if (line.find("#CHROM") == 0)
			{
				headerLines.emplace_back(line);
			}
			else
			{
				headerColumns = line;
			}
		}
		std::string columnLine = setSamplePtrs(line, bamSamplePtrs);
		headerLines.emplace_back(columnLine);
		this->m_vcf_writer->writeHeader(headerLines);
		this->m_vcf_writer->setSamples(columnLine, this->m_sample_ptrs_map);
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
		std::vector< std::string > standardColumnNames = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
		for (auto column : columns)
		{
			std::string upperColumn;
			std::transform(column.begin(), column.end(), std::back_inserter(upperColumn), ::toupper);
			if (std::find(standardColumnNames.begin(), standardColumnNames.end(), upperColumn) == standardColumnNames.end())
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

	/*
	std::string VCFWriter::getVCFColumns(const std::string& line, std::vector< Sample::SharedPtr >& samplePtrs)
	{
		std::string newLine = "";
		std::vector< std::string > sampleNames;
		std::vector< std::string > columns;
		split(line, '\t', columns);
		std::vector< std::string > standardColumnNames = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"};
		for (auto column : columns)
		{
			std::string upperColumn;
			std::transform(column.begin(), column.end(), std::back_inserter(upperColumn), ::toupper);
			if (std::find(standardColumnNames.begin(), standardColumnNames.end(), upperColumn) == standardColumnNames.end())
			{
				sampleNames.emplace_back(column);
			}
			newLine += column + "\t";
		}

		std::vector< Sample::SharedPtr > newSamplePtrs;
		for (auto sampleName : sampleNames)
		{
			for (auto samplePtr : samplePtrs)
			{
				if (sampleName == samplePtr->getName())
				{
					newSamplePtrs.emplace_back(samplePtr);
				}
				else
				{
					auto tmpSamplePtr = std::make_shared< Sample >(sampleName, "", "");
					newSamplePtrs.emplace_back(tmpSamplePtr);
				}
			}
		}
		newLine.pop_back();
		setSamplePtrs(newSamplePtrs);
	}

	void VCFWriter::setSamplePtrs(std::vector< Sample::SharedPtr >& samplePtrs)
	{
		this->m_sample_ptrs_map.clear();
		for (auto samplePtr : samplePtrs)
		{
			if (this->m_sample_ptrs_map.find(samplePtr->getName()) != this->m_sample_ptrs_map.end())
			{
				std::cout << "duplicate sample name: " << samplePtr->getName() << std::endl;
			}
			this->m_sample_ptrs_map.emplace(samplePtr->getName(), samplePtr);
		}
	}

	Sample::SharedPtr VCFWriter::getSamplePtr(const std::string& sampleName)
	{
		auto iter = this->m_sample_ptrs_map.find(sampleName);
		if (iter != this->m_sample_ptrs_map.end())
		{
			return iter->second;
		}
		else
		{
			return nullptr;
		}
	}
	*/
}
