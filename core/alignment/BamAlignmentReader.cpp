#include "BamAlignmentReader.h"
#include "BamAlignment.h"
#include "AlignmentList.h"
#include "core/sample/Sample.h"

#include <unordered_set>

namespace graphite
{
	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath) :
		m_bam_path(bamPath),
		m_is_open(false),
		m_alignment_reader_manager_ptr(nullptr)
	{
	}

	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath, AlignmentReaderManager< BamAlignmentReader >* alignmentReaderManagerPtr) :
	    m_bam_path(bamPath),
		m_is_open(false),
		m_alignment_reader_manager_ptr(alignmentReaderManagerPtr)
	{
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
	}

	void BamAlignmentReader::open()
	{
		std::lock_guard< std::mutex > l(m_lock);
		if (m_is_open) { return; }
		m_is_open = true;
		if (!this->m_bam_reader.Open(this->m_bam_path))
		{
			throw "Unable to open bam file";
		}
		this->m_bam_reader.LocateIndex();
	}

	void BamAlignmentReader::close()
	{
		std::lock_guard< std::mutex > l(m_lock);
		if (!m_is_open) { return; }
		m_is_open = false;
		this->m_bam_reader.Close();
	}

	std::vector< IAlignment::SharedPtr > BamAlignmentReader::loadAlignmentsInRegion(Region::SharedPtr regionPtr, SampleManager::SharedPtr sampleManagerPtr, bool excludeDuplicateReads)
	{
		if (!m_is_open)
		{
			std::cout << "Bam file not opened" << std::endl;
			exit(0);
		}
		std::vector< IAlignment::SharedPtr > alignmentPtrs;

		int refID = this->m_bam_reader.GetReferenceID(regionPtr->getReferenceID());
		// add 1 to the start and end positions because this is 0 based
		this->m_bam_reader.SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());


		size_t counter = 0;

		uint32_t count = 0;
		BamTools::BamAlignment bamAlignment;
		while(this->m_bam_reader.GetNextAlignment(bamAlignment))
		{
            if (bamAlignment.IsDuplicate() && excludeDuplicateReads) { continue; }
			std::string sampleName;
			bamAlignment.GetTag("RG", sampleName);

			Sample::SharedPtr samplePtr = sampleManagerPtr->getSamplePtr(sampleName);
			if (samplePtr == nullptr)
			{
				throw "There was an error in the sample name for: " + sampleName;
			}
			alignmentPtrs.push_back(std::make_shared< BamAlignment >(bamAlignment, samplePtr));
		}
		if (m_alignment_reader_manager_ptr != nullptr)
		{
			m_alignment_reader_manager_ptr->checkinReader(this->shared_from_this());
		}
		return alignmentPtrs;
	}

	std::vector< Sample::SharedPtr > BamAlignmentReader::GetBamReaderSamples(const std::string& bamPath)
	{
		std::vector< Sample::SharedPtr > samplePtrs;
		BamTools::BamReader bamReader;
		if (!bamReader.Open(bamPath))
		{
			throw "Unable to open bam file";
		}
		auto readGroups = bamReader.GetHeader().ReadGroups;
		auto iter = readGroups.Begin();
		for (; iter != readGroups.End(); ++iter)
		{
			auto samplePtr = std::make_shared< Sample >((*iter).Sample, (*iter).ID, bamPath);
			samplePtrs.emplace_back(samplePtr);
		}
		bamReader.Close();
		return samplePtrs;
	}

	position BamAlignmentReader::GetLastPositionInBam(const std::string& bamPath, Region::SharedPtr regionPtr)
	{
		BamTools::BamReader bamReader;
		if (!bamReader.Open(bamPath))
		{
			throw "Unable to open bam file";
		}

		bamReader.LocateIndex();
		int refID = bamReader.GetReferenceID(regionPtr->getReferenceID());
		auto referenceData = bamReader.GetReferenceData();
		bamReader.Close();
		return referenceData[refID].RefLength;
	}

	uint32_t BamAlignmentReader::GetReadLength(const std::string& bamPath)
	{
		uint32_t bamReadLength = 300;
		BamTools::BamReader bamReader;
		if (!bamReader.Open(bamPath))
		{
			throw "Unable to open bam file";
		}
		BamTools::BamAlignment bamAlignment;
		while(bamReader.GetNextAlignment(bamAlignment))
		{
			if (bamAlignment.IsPrimaryAlignment())
			{
				bamReadLength = bamAlignment.QueryBases.size();
				break;
			}
		}
		bamReader.Close();
		return bamReadLength;
	}
}
