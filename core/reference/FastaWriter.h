#ifndef GRAPHITE_FASTAWRITER_H
#define GRAPHITE_FASTAWRITER_H

#include <boost/noncopyable.hpp>
#include <string>
#include <memory>

namespace graphite
{
	class FastaWriter : boost::noncopyable
	{
	public:
		typedef std::shared_ptr< FastaWriter > SharedPtr;
		FastaWriter(const std::string& header, const std::string& sequence);
		~FastaWriter();

		void write(std::ostream& out);

	private:
		std::string m_header;
		std::string m_sequence;
	};
}

#endif //GRAPHITE_FASTAWRITER_H
