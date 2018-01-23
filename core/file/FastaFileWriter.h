#ifndef GRAPHITE_FASTAFILEWRITER_H
#define GRAPHITE_FASTAFILEWRITER_H

#include "ASCIIFileWriter.h"

#include <unordered_set>
#include <mutex>

namespace graphite
{
	class FastaFileWriter : public ASCIIFileWriter
	{
	public:
		typedef std::shared_ptr< FastaFileWriter > SharedPtr;
		FastaFileWriter(const std::string& path);
		virtual ~FastaFileWriter();

		void writeSequence(const std::string& description, const std::string& sequence);

	private:

		std::mutex m_write_lock;
		std::mutex m_descriptions_lock;
		std::unordered_set< std::string > m_descriptions; // this is a lookup table so we can check which desciptions we have written to file
	};
}

#endif //GRAPHITE_FASTAFILEWRITER_H
