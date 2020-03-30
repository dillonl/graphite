#include "AlignmentReader.h"

#include "core/util/Utility.h"

#include "cram/cram.h"
#include "cram/cram_io.h"

namespace graphite
{
	// Example from here: https://www.biostars.org/p/151053/
	AlignmentReader::AlignmentReader(const std::string& filename, const std::string& refPath) : m_path(filename), m_overwrite_sample("")
	{
		m_in = hts_open(m_path.c_str(), "r");
		m_header = sam_hdr_read(m_in);
		hts_set_fai_filename(m_in, refPath.c_str());
		this->m_idx = sam_index_load(m_in, this->m_path.c_str());
		this->init();
	}

	AlignmentReader::~AlignmentReader()
	{
		bam_hdr_destroy(m_header);
		hts_close(m_in);
	}

	void AlignmentReader::init()
	{
		this->m_sample_ptrs.clear();
		this->m_available_regions.clear();
		std::string headerText = std::string(m_header->text);
		std::vector< std::string > readGroups;
		std::vector< std::string > lines;
		split(headerText, '\n', lines);
		for (auto line : lines)
		{
			if (line.size() >= 3 && line.substr(0, 3).compare("@RG") == 0)
			{
				std::vector< std::string > lineSplit;
				split(line, '\t', lineSplit);
				std::string sm;
				std::string id;
				for (std::string lineComponents : lineSplit)
				{
					std::vector< std::string > rgComps;
					split(lineComponents, ':', rgComps);
					if (rgComps.size() > 1 && rgComps[0].compare("SM") == 0)
					{
						sm = rgComps[1];
					}
					if (rgComps.size() > 1 && rgComps[0].compare("ID") == 0)
					{
						id = rgComps[1];
					}
				}
				if (!sm.empty() && !id.empty())
				{
					auto samplePtr = std::make_shared< Sample >(sm, id, m_path);
					this->m_sample_ptrs.emplace(id, samplePtr);
				}
			}
			else if (line.size() >= 3 && line.substr(0, 3).compare("@SQ") == 0)
			{
				std::vector< std::string > lineSplit;
				split(line, '\t', lineSplit);
				for (std::string lineComponents : lineSplit)
				{
					std::vector< std::string > rgComps;
					split(lineComponents, ':', rgComps);
					if (rgComps.size() > 1 && rgComps[0].compare("SN") == 0)
					{
						this->m_available_regions.emplace_back(rgComps[1]);
					}
				}
			}
		}
		setReadLength();
	}

	void AlignmentReader::setReadLength()
	{
		this->m_read_length = 0;

		bam1_t* alignmentPtr = bam_init1();
		for (auto region : this->m_available_regions)
		{
			hts_itr_t* iter = sam_itr_querys(m_idx, m_header, region.c_str());
			auto isSet = sam_itr_next(m_in, iter, alignmentPtr);
			if (isSet > 0)
			{
				this->m_read_length = alignmentPtr->core.l_qseq;
				break;
			}
		}
		if (this->m_read_length == 0)
		{
			std::cout << "There was a problem reading the sample file: " << m_path << std::endl;
			exit(0);
		}
		bam_destroy1(alignmentPtr);
	}

	void AlignmentReader::fetchAlignmentPtrsInRegion(std::vector< std::shared_ptr< Alignment > >& alignmentPtrs, Region::SharedPtr regionPtr, bool unmappedOnly, bool includeDuplicateReads, int32_t mappingQuality)
	{
		bam1_t* htsAlignmentPtr = bam_init1();
		hts_itr_t* iter = sam_itr_querys(this->m_idx, m_header, regionPtr->getRegionString().c_str());

		while ( sam_itr_next(m_in, iter, htsAlignmentPtr) >= 0)
		{
			std::string readGroup = "";
			if (this->m_overwrite_sample.size() > 0)
			{
				auto iter = this->m_sample_ptrs.find(this->m_overwrite_sample);
				readGroup = iter->first;
			}
			else
			{
				std::string readGroup = std::string((char*)bam_aux_get(htsAlignmentPtr, "RG"));
				readGroup = readGroup.substr(1); // the last char is garbage
			}
			position pos = htsAlignmentPtr->core.pos;
			char* name  = bam_get_qname(htsAlignmentPtr);
			bool firstMate = (htsAlignmentPtr->core.flag & BAM_FREAD1);
			bool isMapped = !(htsAlignmentPtr->core.flag & BAM_FUNMAP);
			bool forwardStrand = !bam_is_rev(htsAlignmentPtr);
			bool duplicate = (htsAlignmentPtr->core.flag & BAM_FDUP);
			uint16_t mapQuality = *bam_get_qual(htsAlignmentPtr);
			std::string readName(name);
			auto iter = m_sample_ptrs.find(readGroup);
			if (iter == m_sample_ptrs.end())
			{
				continue;
			}

			auto alignmentPtr = std::make_shared< Alignment >((char*)bam_get_seq(htsAlignmentPtr), htsAlignmentPtr->core.l_qseq, readName, forwardStrand, firstMate, mapQuality, iter->second);
			alignmentPtrs.emplace_back(alignmentPtr);
		}

		sam_itr_destroy(iter);
		bam_destroy1(htsAlignmentPtr);

	}

	void AlignmentReader::overwriteSample(Sample::SharedPtr samplePtr)
	{
		this->m_sample_ptrs.empty();
		this->m_sample_ptrs.emplace(samplePtr->getName(), samplePtr);
		this->m_overwrite_sample = samplePtr->getName();
	}

}
