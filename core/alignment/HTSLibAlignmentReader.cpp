#include "HTSLibAlignmentReader.h"
#include "HTSLibAlignment.h"
#include "SampleManager.hpp"
#include "core/util/Utility.h"

#include <string.h>

namespace graphite
{

	HTSLibAlignmentReader::HTSLibAlignmentReader(const std::string& filePath) :
		m_file_path(filePath),
		m_is_open(false)
	{
	}
	HTSLibAlignmentReader::HTSLibAlignmentReader(const std::string& filePath, AlignmentReaderManager< HTSLibAlignmentReader >* alignmentReaderManagerPtr) :
		m_file_path(filePath),
		m_is_open(false),
		m_alignment_reader_manager_ptr(alignmentReaderManagerPtr)
	{
	}

	HTSLibAlignmentReader::~HTSLibAlignmentReader()
	{
		if (!m_is_open) { close(); }
	}

	void HTSLibAlignmentReader::open()
	{
		std::lock_guard< std::mutex > l(m_lock);
		if (m_is_open) { return; }
		m_is_open = true;
		m_bam_alignment = bam_init1();
		m_fp = sam_open(m_file_path.c_str(), "r");
		m_bam_header = sam_hdr_read(m_fp);
		m_bam_index= bam_index_load(m_file_path.c_str());
	}

	void HTSLibAlignmentReader::close()
	{
		std::lock_guard< std::mutex > l(m_lock);
		if (!m_is_open) { return; }
		m_is_open = false;
		bam_destroy1(m_bam_alignment);
		bam_hdr_destroy(m_bam_header);
		hts_idx_destroy(m_bam_index);
		sam_close(m_fp);
	}

	std::vector< Sample::SharedPtr > HTSLibAlignmentReader::GetBamReaderSamples(const std::string& bamPath)
	{
		// a quick note to future people that are planning on looking at this function:
		// this is a horrible function that is just chores. There are no cool
		// algorithmic fireworks here, just brute-force string splitting.
		std::vector< Sample::SharedPtr > samplePtrs;
		auto readerPtr = std::make_shared< HTSLibAlignmentReader >(bamPath);
		readerPtr->open();

		std::vector< std::string > headerLines;
		split(readerPtr->m_bam_header->text, '\n', headerLines);
		for (auto headerLine : headerLines)
		{
			if (headerLine.compare(0, 3, "@RG") == 0)
			{
				std::vector< std::string > readGroupLine;
				split(headerLine, '\t', readGroupLine);
				std::string readGroupID = "";
				std::string readGroupSampleName = "";
				for (auto readGroupInfo : readGroupLine)
				{
					std::vector< std::string > infoSplit;
					split(readGroupInfo, ':', infoSplit);
					if (infoSplit.size() < 2) { continue; }
					if (infoSplit[0].compare(0, 2, "ID") == 0)
					{
						readGroupID = infoSplit[1];
					}
					else if (infoSplit[0].compare(0, 2, "SM") == 0)
					{
						readGroupSampleName = infoSplit[1];
					}
				}
				if (readGroupID.size() > 0 && readGroupSampleName.size() > 0)
				{
					auto samplePtr = std::make_shared< Sample >(readGroupSampleName, readGroupID, bamPath);
					samplePtrs.emplace_back(samplePtr);
				}
			}
		}

		readerPtr->close();
		return samplePtrs;
	}

	std::vector< IAlignment::SharedPtr > HTSLibAlignmentReader::loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads)
	{
		std::lock_guard< std::mutex > l(m_lock);
		if (!m_is_open)
		{
			std::cout << "file must be open before you read the region" << std::endl;
			exit(0);
		}
		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		// std::vector< Roaring > buckets(m_lsh_ptr->getNumBuckets());
		hts_itr_t* iter = sam_itr_querys(m_bam_index, m_bam_header, regionPtr->getRegionString().c_str());
		char* qseq = (char*)malloc(1024);
		while (sam_itr_next(m_fp, iter, m_bam_alignment) >= 0)
		{
			if (m_bam_alignment->core.flag & 0x900 || (excludeDuplicateReads && m_bam_alignment->core.flag & 0x400)) { continue; }
			uint32_t len = m_bam_alignment->core.l_qseq; //length of the read.

			Sample::SharedPtr samplePtr = SampleManager::Instance()->getSamplePtr(bam_aux2Z(bam_aux_get(m_bam_alignment, "RG")));
			if (samplePtr == nullptr)
			{
				throw "There was an error in the sample name for: " + std::string(bam_aux2Z(bam_aux_get(m_bam_alignment, "RG")));
			};
			auto alignmentPtr = std::make_shared< HTSLibAlignment >();
			alignmentPtr->setPosition(m_bam_alignment->core.pos);
			alignmentPtr->setSample(samplePtr);
			alignmentPtr->setIsReverseStrand(bam_is_rev(m_bam_alignment));
			alignmentPtr->setName(bam_get_qname(m_bam_alignment), (m_bam_alignment->core.flag & 0x0040));
			alignmentPtr->setDuplicate(m_bam_alignment->core.flag & 0x00400);
			alignmentPtr->setFilePosition(bgzf_tell(m_fp->fp.bgzf));
			alignmentPtrs.emplace_back(alignmentPtr);

			uint8_t* q = bam_get_seq(m_bam_alignment); //quality string
			for(int i=0; i < len; ++i)
			{
				qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
			}
			alignmentPtr->setSequence(qseq, len);

			/*
			AlignmentSignature sig(0);
			if (AlignmentParser::GetAlignmentSignature(qseq, len, sig))
			{
				m_lsh_ptr->hashSignature(sig, bgzf_tell(m_fp->fp.bgzf));
				// for (auto idx : bucketIdxs)
				// {
					// buckets[idx].emplace_back(bgzf_tell(m_fp->fp.bgzf));
					// buckets[idx].add(bgzf_tell(m_fp->fp.bgzf));
				// }
				++m_alignments_processed;
			}
			*/
		}
		free(qseq);
		hts_itr_destroy(iter);

		// if (auto ptr = m_alignment_reader_manager_ptr.lock())
		// {
			// ptr->checkinReader(this->shared_from_this());
		// }
		m_alignment_reader_manager_ptr->checkinReader(this->shared_from_this());

		return alignmentPtrs;
	}

	void HTSLibAlignmentReader::setAlignmentSequence(uint64_t filePosition, std::string& sequence)
	{
		std::lock_guard< std::mutex > l(m_lock);
		// hts_itr_t* iter = sam_itr_querys(m_bam_index, m_bam_header, regionPtr->getRegionString().c_str());

		if (bgzf_seek(m_fp->fp.bgzf, filePosition, SEEK_SET) == 0 && bam_read1(m_fp->fp.bgzf, m_bam_alignment)) // sam_itr_next(m_fp, iter, m_bam_alignment) >= 0)
		{
			auto len = m_bam_alignment->core.l_qseq;
			sequence.reserve(len + 1);
			uint8_t* q = bam_get_seq(m_bam_alignment); //quality string
			for(int i=0; i < len; ++i)
			{

				 sequence[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
			}
			sequence[len + 1] = '\0';
		}
	}

	position HTSLibAlignmentReader::getRegionLastPosition(Region::SharedPtr regionPtr)
	{
		std::vector< std::string > headerLines;
		split(m_bam_header->text, '\n', headerLines);
		for (auto headerLine : headerLines)
		{
			if (headerLine.compare(0, 3, "@SQ") != 0) { continue; }
			if (headerLine.find("SN:" + regionPtr->getReferenceID()) != std::string::npos)
			{
				auto lengthStart = headerLine.find("LN:") + 3;
				if (lengthStart != std::string::npos)
				{
					auto lengthEnd = headerLine.find("\t", lengthStart);
					if (lengthStart != std::string::npos)
					{
						return stoul(headerLine.substr(lengthStart, lengthEnd - lengthStart));
					}
				}
			}
		}
		return 300000; // just return a large value since we can't find a good match
	}
}
