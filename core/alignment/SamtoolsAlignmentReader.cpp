#include "SamtoolsAlignmentReader.h"
#include "BamAlignment.h"
#include "AlignmentList.h"
#include "SampleManager.hpp"

#include <unordered_set>

#include <boost/algorithm/string.hpp>

#include "bam.h"
#include "sam.h"

namespace graphite
{
	std::unordered_map< std::string, position > SamtoolsAlignmentReader::s_region_last_positions;
	std::unordered_map< std::string, position > SamtoolsAlignmentReader::s_region_map;
	std::mutex SamtoolsAlignmentReader::s_region_last_positions_lock;
	std::mutex SamtoolsAlignmentReader::s_region_ids_lock;

	SamtoolsAlignmentReader::SamtoolsAlignmentReader(const std::string& path) :
		m_path(path)
	{
		InitReader(path); // statically initializes last positions and region information
	}

	SamtoolsAlignmentReader::~SamtoolsAlignmentReader()
	{
	}

	int readAlignment(const bam1_t* b, void* data)
	{
		std::string readGroup = std::string((char*)bam_aux_get(b, "RG"));
		position pos = b->core.pos;
		char* name  = bam1_qname(b);
		bool firstMate = (b->core.flag & BAM_FREAD1);
		bool isMapped = !(b->core.flag & BAM_FUNMAP);
		bool isReverseStrand = bam_is_rev(b);
		bool duplicate = (b->core.flag & BAM_FDUP);
		uint16_t mapQuality = bam_get_qual(b);

		Sample::SharedPtr samplePtr = SampleManager::Instance()->getSamplePtr(readGroup.substr(1)); // I'm not sure why I have to drop the first char
		if (samplePtr == nullptr)
		{
			throw "There was an error in the sample name for: " + std::string(readGroup);
		}

		int n=0;
		char* qseq = (char*)malloc(b->core.l_qseq+1);
		char* s   = bam1_seq(b);
		for(n=0;n<(b->core.l_qseq);n++)
		{
			char v = bam1_seqi(s,n);
			qseq[n] = bam_nt16_rev_table[v];
		}
		qseq[n] = 0;

		auto bamAlignmentPtr = std::make_shared< BamAlignment >(pos, firstMate, isMapped, isReverseStrand, duplicate, mapQuality, name, samplePtr);
		bamAlignmentPtr->setSequence(qseq, b->core.l_qseq);

		std::vector< IAlignment::SharedPtr >* als = (std::vector< IAlignment::SharedPtr >*)data;
		als->emplace_back(bamAlignmentPtr);
		return 0;
	}

	int readAlignmentSequences(const bam1_t* b, void* data)
	{
		std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > nameAlignmentPtrsMap = (*(std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > >*)data);
		bool firstMate = (b->core.flag & BAM_FREAD1);
		std::string name  = std::string(bam1_qname(b) + std::to_string(firstMate));

		auto iter = nameAlignmentPtrsMap->find(name);
		if (iter == nameAlignmentPtrsMap->end())
		{
			return 0;
		}
		IAlignment::SharedPtr alignmentPtr = iter->second;
		int n=0;
		char* qseq = (char*)malloc(b->core.l_qseq+1);
		char* s   = bam1_seq(b);
		for(n=0;n<(b->core.l_qseq);n++)
		{
			char v = bam1_seqi(s,n);
			qseq[n] = bam_nt16_rev_table[v];
		}
		qseq[n] = 0;
		alignmentPtr->setSequence(qseq, b->core.l_qseq);
		// position pos = b->core.pos;
		// std::cout << pos << " " << qseq << std::endl;
		return 0;
	}

	typedef struct {
		int beg, end;
		samfile_t *in;
	} tmpstruct_t;

	// callback for bam_plbuf_init()
	static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
	{
		tmpstruct_t *tmp = (tmpstruct_t*)data;
		// if ((int)pos >= tmp->beg && (int)pos < tmp->end)
			// printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
		return 0;
	}

	std::vector< IAlignment::SharedPtr > SamtoolsAlignmentReader::loadAlignmentsInRegion(Region::SharedPtr regionPtr, bool excludeDuplicateReads)
	{
		std::vector< IAlignment::SharedPtr > alignmentPtrs;

		int ref;
		int beg = 0;
		int end  = 0x7fffffff;

		static std::mutex l;
		std::lock_guard< std::mutex > lo(l);
		auto fp = samopen(m_path.c_str(), "r", 0);

		auto idx = bam_index_load(m_path.c_str());

		bam_plbuf_t* buf;
		void* data = &alignmentPtrs;
		bam_parse_region(fp->header, regionPtr->getRegionString().c_str(), &ref, &beg, &end); // parse the region
		buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
		bam_fetch(fp->x.bam, idx, ref, beg, end, data, readAlignment);

		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
		samclose(fp);

		return alignmentPtrs;
	}

	void SamtoolsAlignmentReader::loadAlignmentSequencesInRegion(Region::SharedPtr regionPtr, std::shared_ptr< std::unordered_map< std::string, IAlignment::SharedPtr > > nameAlignmentPtrsMap)
	{
		int ref;
		tmpstruct_t tmp;
		tmp.beg = 0;
		tmp.end = 0x7fffffff;
		tmp.in = samopen(m_path.c_str(), "r", 0);
		auto idx = bam_index_load(m_path.c_str());

		bam_plbuf_t* buf;
		void* data = &nameAlignmentPtrsMap;
		bam_parse_region(tmp.in->header, regionPtr->getRegionString().c_str(), &ref, &tmp.beg, &tmp.end); // parse the region
		buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
		bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, data, readAlignmentSequences);

		bam_plbuf_push(0, buf); // finalize pileup
		bam_index_destroy(idx);
		bam_plbuf_destroy(buf);
	}

	std::vector< Sample::SharedPtr > SamtoolsAlignmentReader::GetBamReaderSamples(const std::string& path)
	{
		std::vector< Sample::SharedPtr > samplePtrs;
		auto fp = bam_open(path.c_str(), "r");
		bam_hdr_t* header = bam_header_read(fp);
		auto idx = bam_index_load(path.c_str());

		std::string headerText = std::string(header->text);
		std::vector< std::string > readGroups;
		std::vector< std::string > lines;
		boost::split(lines, headerText, boost::is_any_of("\n"));
		for (auto line : lines)
		{
			if (line.size() >= 3 && line.substr(0, 3).compare("@RG") == 0)
			{
				std::vector< std::string > lineSplit;
				boost::split(lineSplit, line, boost::is_any_of("\t"));
				std::string sm;
				std::string id;
				for (std::string lineComponents : lineSplit)
				{
					std::vector< std::string > rgComps;
					boost::split(rgComps, lineComponents, boost::is_any_of(":"));
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
					auto samplePtr = std::make_shared< Sample >(sm, id, path);
					samplePtrs.emplace_back(samplePtr);
				}
			}
		}

		bam_header_destroy(header);
		bam_index_destroy(idx);
		bam_close(fp);
		return samplePtrs;
	}

	position SamtoolsAlignmentReader::GetLastPositionInBam(const std::string& bamPath, Region::SharedPtr regionPtr)
	{
		std::lock_guard< std::mutex > lock(s_region_last_positions_lock);
		return s_region_last_positions[regionPtr->getReferenceID()];
	}

	void SamtoolsAlignmentReader::InitReader(const std::string& path)
	{
		std::lock_guard< std::mutex > lock(s_region_last_positions_lock);
		if (s_region_last_positions.size() > 0) { return; }
		auto fp = bam_open(path.c_str(), "r");
		/* bgzf_set_cache_size(fp, 8 * 1024 *1024); */
		bam_hdr_t* header = bam_header_read(fp);
		auto idx = bam_index_load(path.c_str());

		for (auto i = 0; i < header->n_targets; ++i)
		{
			s_region_last_positions.emplace(std::string(header->target_name[i]), header->target_len[i]);
			s_region_map.emplace(std::string(header->target_name[i]), i);
		}

		bam_header_destroy(header);
		bam_index_destroy(idx);
		bam_close(fp);
	}
}
