#include "BamFileWriter.h"

#include "api/BamReader.h"
// #include "api/Sort.h"
#include "api/algorithms/Sort.h"

#include <cstdio>
#include <unordered_set>

namespace graphite
{
	BamFileWriter::BamFileWriter(const std::string& path, BamTools::SamHeader header, BamTools::RefVector& refVec)
		: m_path(path), m_header(header), m_ref_vec(refVec), m_open(true)
	{
        m_tmp_path = m_path + ".tmp";
		m_tmp_bam_writer.Open(m_tmp_path, m_header, m_ref_vec);
		m_tmp_bam_writer.SetCompressionMode(BamTools::BamWriter::CompressionMode::Uncompressed);
	}

	BamFileWriter::~BamFileWriter()
	{
		close();
	}

	void BamFileWriter::writeAlignment(std::shared_ptr< BamTools::BamAlignment > bamAlignmentPtr, const std::string& ref)
	{
		std::lock_guard< std::mutex > l(m_writer_lock);
		if (m_ref_map.find(ref) == m_ref_map.end())
		{
			m_ref_id_map.emplace(ref, m_ref.size());  // add the idx for the region for fast lookup
			m_ref_map.emplace(ref, 1);
			m_ref.emplace_back(ref);
		}
		else
		{
			m_ref_map[ref] += 1;
		}
		bamAlignmentPtr->RefID = m_ref_id_map[ref] + m_header.Sequences.Size();
		m_tmp_bam_writer.SaveAlignment(*bamAlignmentPtr);
	}

	void printRefVec(std::vector< BamTools::CigarOp >& refVec)
	{
		for (auto rv : refVec)
		{
			std::cout << rv.Type << std::to_string(rv.Length);
		}
		std::cout << std::endl;
	}

	void BamFileWriter::close()
	{
		std::lock_guard< std::mutex > l(m_file_open_lock);
		if (m_open)
		{
			m_open = false;
			for (auto i = 0; i <  m_ref_vec.size(); ++i)
			{
				m_ref_vec[i].RefLength = 0;
			}
			std::vector< BamTools::SamSequence > sequences;
			for (auto iter = m_header.Sequences.Begin(); iter != m_header.Sequences.End(); ++iter)
			{
				BamTools::SamSequence seq(iter->Name, 0);
				sequences.emplace_back(seq);
			}
			m_header.Sequences.Clear();
			for (auto seq : sequences)
			{
				m_header.Sequences.Add(seq);
			}
			for (auto ref : m_ref)
			{
				BamTools::SamSequence s(ref,m_ref_map[ref]);
				m_header.Sequences.Add(s);
				BamTools::RefData rd(ref, m_ref_map[ref]);
				m_ref_vec.emplace_back(rd);
			}
            m_tmp_bam_writer.Close();

			{
				BamTools::BamWriter bamWriter;
				bamWriter.Open(m_path, m_header, m_ref_vec);

				{
					BamTools::BamReader bamReader;
					if (!bamReader.Open(m_tmp_path))
					{
						throw "Unable to open bam file";
					}

					std::vector< BamTools::BamAlignment > alignments;
					BamTools::BamAlignment bamAlignment;
					while (bamReader.GetNextAlignment(bamAlignment))
					{
						alignments.emplace_back(bamAlignment);
					}

					std::sort(alignments.begin(), alignments.end(), BamTools::Algorithms::Sort::ByPosition());
					int count = 0;
					for (auto bamAlignment : alignments) // while true is never a good answer, maybe I'll rethink this in the future
					{
						bamWriter.SaveAlignment(bamAlignment);
					}

					bamReader.Close();

				}

				bamWriter.Close();
			}
			{
				BamTools::BamReader bamReader;
				if (!bamReader.Open(m_path))
				{
					throw "Unable to open bam file";
				}
				if (!bamReader.CreateIndex())
				{
					std::cout << "failed to create index" << std::endl;
				}
				bamReader.Close();
                remove(m_tmp_path.c_str());
            }
		}
	}
}
