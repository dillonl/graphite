#include "BamReader.h"

#include <unordered_map>

namespace graphite
{
	BamReader::BamReader(const std::string& filename) : m_bam_path(filename),
														m_overwrite_sample(false)
	{
		this->m_bam_reader = std::make_shared< BamTools::BamReader >();
		if (!this->m_bam_reader->Open(filename))
		{
			throw "Unable to open bam file";
		}
		this->m_bam_reader->LocateIndex();
		initializeSamplePtrs();
	}

	BamReader::~BamReader()
	{
		this->m_bam_reader->Close();
	}

	void BamReader::initializeSamplePtrs()
	{
		auto readGroups = this->m_bam_reader->GetHeader().ReadGroups;
		auto iter = readGroups.Begin();
		for (; iter != readGroups.End(); ++iter)
		{
			auto samplePtr = std::make_shared< Sample >((*iter).Sample, (*iter).ID, this->m_bam_path);
			this->m_sample_ptrs.emplace(samplePtr);
		}
	}

	void BamReader::overwriteSample(Sample::SharedPtr samplePtr)
	{
		this->m_sample_ptrs.empty();
		this->m_sample_ptrs.emplace(samplePtr);
		this->m_overwrite_sample = true;
	}

	std::unordered_set< Sample::SharedPtr > BamReader::getSamplePtrs()
	{
		return this->m_sample_ptrs;
	}

	/*
	 * Concat new bamAlignments to the passed in list, bamAlignmentPtrs.
	 */
	void BamReader::fetchBamAlignmentPtrsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs,  Region::SharedPtr regionPtr, bool unmappedOnly, bool includeDuplicateReads, int32_t mappingQuality)
	{
		int refID = this->m_bam_reader->GetReferenceID(regionPtr->getReferenceID());
		this->m_bam_reader->SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());
		BamTools::BamAlignment* bamtoolsAlignmentPtr = new BamTools::BamAlignment();
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();
		while (this->m_bam_reader->GetNextAlignment(*bamtoolsAlignmentPtr))
		{
			std::string sampleName;
			bamtoolsAlignmentPtr->GetTag("RG", sampleName);
			bool isInRegion = (startPosition < bamtoolsAlignmentPtr->Position && (bamtoolsAlignmentPtr->Position + bamtoolsAlignmentPtr->Length) < endPosition);
			bool shouldFilter = (bamtoolsAlignmentPtr->MapQuality <= mappingQuality); // remove mapping quality less than param value (default is -1 aka no filter)
			if ((bamtoolsAlignmentPtr->IsDuplicate() && !includeDuplicateReads) ||
				(unmappedOnly && bamtoolsAlignmentPtr->IsMapped()) || !isInRegion || shouldFilter)
			{
				delete bamtoolsAlignmentPtr; // delete the ptr if the ptr isn't added to the alignmentPtrs
			}
			else
			{
				std::shared_ptr< BamTools::BamAlignment > bamAlignmentSharedPtr(bamtoolsAlignmentPtr, [](BamTools::BamAlignment* btPtr)
																				{
																					delete btPtr;
																				});
				bamAlignmentPtrs.emplace_back(bamAlignmentSharedPtr);
			}
			bamtoolsAlignmentPtr = new BamTools::BamAlignment();
		}
		delete bamtoolsAlignmentPtr; // the last bamtoolsAlignmentPtr is always when GetNextAlignment returned false so we can delete it
	}

	uint32_t BamReader::getReadLength()
	{
		BamTools::BamAlignment* bamtoolsAlignmentPtr = new BamTools::BamAlignment();
		this->m_bam_reader->GetNextAlignment(*bamtoolsAlignmentPtr);
		uint32_t readLength = bamtoolsAlignmentPtr->Length;
		delete bamtoolsAlignmentPtr;
		return readLength;
	}
}
