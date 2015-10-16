#ifndef GRAPHITE_VCFHEADER_H
#define GRAPHITE_VCFHEADER_H

#include <memory>
#include <vector>
#include <string>

#include <boost/noncopyable.hpp>

namespace graphite
{
	class VCFHeader : private boost::noncopyable
	{
	public:
		typedef std::shared_ptr< VCFHeader > SharedPtr;
		VCFHeader();
		~VCFHeader();

		void addHeaderLine(const std::string& headerLine);
		std::string getHeader();
	private:
		std::vector< std::string > m_header_lines;
	};
}

#endif //GRAPHITE_VCFHEADER_H
