#ifndef GRAPHITE_FASTAFILEWRITER_H
#define GRAPHITE_FASTAFILEWRITER_H

//#include <string>
//#include <memory>
//#include <ostream>
#include <vector>
#include <fstream>

/*
 * Why inherit Noncopyable.hpp?
 * What is the typdef line for?
 */

#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class FastaFileWriter //: Noncopyable
	{
        /*
	public:
		typedef std::shared_ptr< FastaWriter > SharedPtr;
		FastaWriter(const std::string& header, const std::string& sequence);
		~FastaWriter();

		void write(std::ostream& out);

	private:
		std::string m_header;
		std::string m_sequence;
        */

    public:
        //typedef std::shared_ptr< FastaFileWriter > SharedPtr;
        FastaFileWriter (std::vector< std::string >& headers, std::vector< std::string >& sequences);

        ~FastaFileWriter ();

        bool open (std::string fileName);
        void write ();
        bool close ();

    private:
        bool m_opened;
        std::ofstream m_out_stream;
        std::vector< std::string > m_headers;
        std::vector< std::string > m_sequences;
	};
}

#endif //GRAPHITE_FASTAWRITER_H
