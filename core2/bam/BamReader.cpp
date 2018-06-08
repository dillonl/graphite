#include "BamReader.h"

namespace graphite
{
	BamReader::BamReader(const std::string& filename) : m_bam_path(filename)
	{
		this->m_bam_reader = std::make_shared< BamTools::BamReader >();
		if (!this->m_bam_reader->Open(filename))
		{
			throw "Unable to open bam file";
		}
		this->m_bam_reader->LocateIndex();
	}

	BamReader::~BamReader()
	{
		this->m_bam_reader->Close();
	}

	std::vector< Sample::SharedPtr > BamReader::getSamplePtrs()
	{
		std::vector< Sample::SharedPtr > samplePtrs;
		auto readGroups = this->m_bam_reader->GetHeader().ReadGroups;
		auto iter = readGroups.Begin();
		for (; iter != readGroups.End(); ++iter)
		{
			auto samplePtr = std::make_shared< Sample >((*iter).Sample, (*iter).ID, this->m_bam_path);
			samplePtrs.emplace_back(samplePtr);
		}
		return samplePtrs;
	}

	/*
	 * Concat new bamAlignments to the passed in list, bamAlignmentPtrs.
	 */
	void BamReader::fetchBamAlignmentPtrsInRegion(std::vector< std::shared_ptr< BamAlignment > >& bamAlignmentPtrs,  Region::SharedPtr regionPtr, bool unmappedOnly, bool includeDuplicateReads)
	{
		int refID = this->m_bam_reader->GetReferenceID(regionPtr->getReferenceID());
		this->m_bam_reader->SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());
		BamTools::BamAlignment* bamtoolsAlignmentPtr = new BamTools::BamAlignment();
		while (this->m_bam_reader->GetNextAlignment(*bamtoolsAlignmentPtr))
		{

			if ((bamtoolsAlignmentPtr->IsDuplicate() && !includeDuplicateReads) ||
				(unmappedOnly && bamtoolsAlignmentPtr->IsMapped()))
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
		return 0;
	}
}
