#include "BamAlignmentReaderPreloadManager.h"
#include "core/util/ThreadPool.hpp"

namespace gwiz
{
	BamAlignmentReaderPreloadManager::BamAlignmentReaderPreloadManager(const std::string& bamPath) :
		m_bam_path(bamPath), m_region(nullptr), m_alignments_ptr(std::make_shared< std::vector< BamAlignment::SharedPtr > >())
	{

		std::cout << this->m_alignments_ptr.get() << std::endl;
		// processBam();
	}

	BamAlignmentReaderPreloadManager::BamAlignmentReaderPreloadManager(const std::string& bamPath, Region::SharedPtr region) :
		m_bam_path(bamPath), m_region(region), m_alignments_ptr(std::make_shared< std::vector< BamAlignment::SharedPtr > >())
	{
		// processBam();
	}

	BamAlignmentReaderPreloadManager::~BamAlignmentReaderPreloadManager()
	{
	}

	IAlignmentReader::SharedPtr BamAlignmentReaderPreloadManager::generateAlignmentReader()
	{
		return std::make_shared< BamAlignmentReaderPreload >(this->m_alignments_ptr);
	}

	/*
	void BamAlignmentReaderPreloadManager::processBam()
	{
		if (!this->m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		if (this->m_region != nullptr)
		{
			this->m_bam_reader.LocateIndex();
			int refID = this->m_bam_reader.GetReferenceID(this->m_region->getReferenceID());
			// add 1 to the start and end positions because this is 0 based
			// this->m_bam_reader.SetRegion(refID, this->m_region->getStartPosition(), refID, this->m_region->getEndPosition());

			// std::vector< std::future< std::vector< std::shared_ptr< BamAlignment > > > > alignmentFutures;
			std::vector< std::shared_ptr< std::future< std::vector< std::shared_ptr< BamAlignment > > > >  > alignmentFuturePtrs;
			position delta = 100000;
			position startPosition = this->m_region->getStartPosition();
			position position = (startPosition <= this->m_region->getEndPosition() - delta) ? startPosition + delta : this->m_region->getEndPosition() - startPosition;
			while (position < this->m_region->getEndPosition())
			{
				uint32_t positionPlusDelta = position + delta;
				uint32_t endPosition = (positionPlusDelta > this->m_region->getEndPosition()) ? this->m_region->getEndPosition() : positionPlusDelta;
				Region::SharedPtr regionPtr = std::make_shared< Region >(this->m_region->getReferenceID() + ":" + std::to_string(position) + "-" + std::to_string(endPosition));
				auto funct = std::bind(&BamAlignmentReaderPreloadManager::getAlignmentsInRegion, this, regionPtr);
				auto functFuture = ThreadPool::Instance()->enqueue(funct);
				alignmentFuturePtrs.emplace_back(functFuture);
				position += delta;
			}
			std::unordered_map< std::string, bool > uniqueAlignments;
			for (auto futureVectorOfAlignmentsPtr : alignmentFuturePtrs)
			{
				futureVectorOfAlignmentsPtr->wait();
				auto vectorOfAlignments = futureVectorOfAlignmentsPtr->get();
				for (auto alignmentPtr : vectorOfAlignments)
				{
					if (uniqueAlignments.find(alignmentPtr->getID()) != uniqueAlignments.end()) { continue; } // if the alignment has been seen then don't add it
					uniqueAlignments[alignmentPtr->getID()] = true;
					this->m_alignments_ptr->emplace_back(alignmentPtr);
				}
			}
			std::sort(this->m_alignments_ptr->begin(), this->m_alignments_ptr->end(),
					  [](const BamAlignment::SharedPtr& a, const BamAlignment::SharedPtr& b)
					  {
						  return a->getPosition() < b->getPosition();
					  });
		}
		else
		{
			BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
			while(this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr))
			{
				// if (counter++ > 1000) { break; }
				this->m_alignments_ptr->push_back(std::make_shared< BamAlignment >(bamAlignmentPtr));
				bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
			}
		}
	}
	*/

	void BamAlignmentReaderPreloadManager::processBam()
	{
		if (!this->m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		if (this->m_region != nullptr)
		{
			this->m_bam_reader.LocateIndex();
			int refID = this->m_bam_reader.GetReferenceID(this->m_region->getReferenceID());
			// add 1 to the start and end positions because this is 0 based
			this->m_bam_reader.SetRegion(refID, this->m_region->getStartPosition(), refID, this->m_region->getEndPosition());
		}
		BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		size_t counter = 0;
		while(this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr))
		{
			// if (counter++ > 1000) { break; }
			this->m_alignments_ptr->push_back(std::make_shared< BamAlignment >(bamAlignmentPtr));
			bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		}
	}

	std::vector< BamAlignment::SharedPtr > BamAlignmentReaderPreloadManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		std::vector< BamAlignment::SharedPtr > alignments;
		BamTools::BamReader bamReader;
		if (!bamReader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		if (regionPtr == nullptr) { throw "In BamAlignmentReaderPreloadManager.cpp -> getAlignmentInRegion must be called with a valid region"; }

		bamReader.LocateIndex();
		int refID = bamReader.GetReferenceID(regionPtr->getReferenceID());
		// add 1 to the start and end positions because this is 0 based
		bamReader.SetRegion(refID, regionPtr->getStartPosition(), refID, regionPtr->getEndPosition());

		BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		while(bamReader.GetNextAlignment(*bamAlignmentPtr))
		{
			alignments.push_back(std::make_shared< BamAlignment >(bamAlignmentPtr));
			bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		}
		return alignments;
	}
}
