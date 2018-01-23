#include "FastaFileWriter.h"

#include "Fasta.h"

namespace graphite
{
	FastaFileWriter::FastaFileWriter(const std::string& path) :
		ASCIIFileWriter(path)
	{
		ASCIIFileWriter::open();
	}

	FastaFileWriter::~FastaFileWriter()
	{
		ASCIIFileWriter::close();
		FastaIndex* fai = new FastaIndex();
		//cerr << "generating fasta index file for " << fastaFileName << endl;
		fai->indexReference(m_file_path);
		std::string fastaIdxPath = m_file_path + fai->indexFileExtension();
		fai->writeIndexFile(fastaIdxPath);
	}

	void FastaFileWriter::writeSequence(const std::string& description, const std::string& sequence)
	{
		{
			std::lock_guard< std::mutex > l(m_descriptions_lock);
			if (m_descriptions.find(description) != m_descriptions.end())
			{
				return; // do not proceed if we have already written this description
			}
			m_descriptions.emplace(description);
		}
		{
			std::lock_guard< std::mutex > l(m_write_lock);
			ASCIIFileWriter::write(">", 1);
			ASCIIFileWriter::write(description.c_str(), description.size());
			ASCIIFileWriter::write("\n", 1);
			ASCIIFileWriter::write(sequence.c_str(), sequence.size());
			ASCIIFileWriter::write("\n", 1);
		}
	}
}
