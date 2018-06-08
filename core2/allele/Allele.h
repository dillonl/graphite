#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core2/util/Noncopyable.hpp"

#include <memory>

namespace graphite
{
	class Allele : private Noncopyable
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
		Allele(const std::string& sequence);
		~Allele();

		std::string getSequence() { return this->m_sequence; }

	private:
		std::string m_sequence;

	};
}

#endif// GRAPHITE_ALLELE_H
