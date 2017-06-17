#ifndef GRAPHITE_FASTAFILEWRITER_H
#define GRAPHITE_FASTAFILEWRITER_H

#include <vector>
#include <fstream>

/*
 * Why inherit Noncopyable.hpp?
 * What is the typdef line for?
 *
 * Add const to constructor
 * Add const to write fxn in cpp file.
 */

#include "IFileWriter.h"

namespace graphite
{
	class FastaFileWriter // : public IFileWriter
	{
    public:
        typedef std::shared_ptr< FastaFileWriter > SharedPtr;
        //FastaFileWriter (const std::string& path);
        FastaFileWriter ();

        ~FastaFileWriter ();

        bool open (std::string fileName);
        //void write (const std::vector< std::string >& headers, const std::vector< std::string >& sequences);
        void write (const std::string& header, const std::string& sequence);
        void close ();
        // std::string getFilePath() overried { return m_file_path; }

    private:
        std::string m_file_path;
        bool m_opened;
        std::ofstream m_out_stream;
	};
}

#endif //GRAPHITE_FASTAWRITER_H
