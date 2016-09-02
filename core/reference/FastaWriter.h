#ifndef GRAPHITE_FASTAWRITER_H
#define GRAPHITE_FASTAWRITER_H

#include <string>
#include <memory>
#include <ostream>

#include "core/util/Noncopyable.hpp"

namespace graphite
{
	class FastaWriter : Noncopyable
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
