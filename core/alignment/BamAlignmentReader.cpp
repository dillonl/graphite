#include "BamAlignmentReader.h"
#include "BamAlignment.h"
#include "AlignmentList.h"
#include "SampleManager.hpp"

#include <unordered_set>

namespace graphite
{
	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath) :
		m_bam_path(bamPath)
	{
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
	}

	std::vector< IAlignment::SharedPtr > BamAlignmentReader::loadAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		// static std::mutex lock;
		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		if (!this->m_bam_reader.Open(this->m_bam_path))
		{
			throw "Unable to open bam file";
		}

		this->m_bam_reader.LocateIndex();
		int refID = this->m_bam_reader.GetReferenceID(regionPtr->getReferenceID());
		// add 1 to the start and end positions because this is 0 based
		this->m_bam_reader.SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());
		// lock.lock();
		// std::cout << "bar region: " << refID << " " << regionPtr->getStartPosition() << " " << regionPtr->getEndPosition();


		// auto bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		size_t counter = 0;

		uint32_t count = 0;
		BamTools::BamAlignment bamAlignment;
		// while(this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr))
		while(this->m_bam_reader.GetNextAlignment(bamAlignment))
		{
			// if (bamAlignment.RefID != refID) { break; }
			std::string sample;
			bamAlignment.GetTag("RG", sample);

			Sample::SharedPtr samplePtr = SampleManager::Instance()->getSamplePtr(sample);
			if (samplePtr == nullptr)
			{
				throw "There was an error in the sample name for: " + sample;
			}
			alignmentPtrs.push_back(std::make_shared< BamAlignment >(bamAlignment, samplePtr));
		}
		this->m_bam_reader.Close();
		// lock.lock();
		// std::cout << " alignments: " << alignmentPtrs.size() << std::endl;
		// lock.unlock();
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
}
