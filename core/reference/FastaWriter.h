#ifndef GWIZ_FASTAWRITER_H
#define GWIZ_FASTAWRITER_H

#include <boost/noncopyable.hpp>
#include <string>

namespace gwiz
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

#endif //GWIZ_FASTAWRITER_H
