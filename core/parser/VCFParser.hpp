#ifndef GRAPHITE_VCFPARSER_H
#define GRAPHITE_VCFPARSER_H

/*
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>

#include "core/sequence/SequenceManager.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace graphite
{

	template<typename Iterator>
	struct VariantParser : boost::spirit::qi::grammar<Iterator, boost::fusion::vector< std::string&, uint32_t&, std::string&, std::string&, std::vector<std::string>&, std::string&, std::string&, std::string& >() >
	{
		VariantParser() : VariantParser::base_type(start)
			{
				using namespace boost;
				using namespace boost::spirit;
				using namespace boost::phoenix;


			}


		boost::spirit::qi::rule<Iterator, void()> tab;
		boost::spirit::qi::rule<Iterator, std::string()> chrom;
		boost::spirit::qi::rule<Iterator, uint32_t()> pos;
		boost::spirit::qi::rule<Iterator, std::string()> id;
		boost::spirit::qi::rule<Iterator, std::string()> ref;
		boost::spirit::qi::rule<Iterator, std::vector<std::string>() > alt;
		boost::spirit::qi::rule<Iterator, std::string()> qual;
		boost::spirit::qi::rule<Iterator, std::string()> filter;
		boost::spirit::qi::rule<Iterator, std::string() > info;

		boost::spirit::qi::rule<Iterator, boost::fusion::vector< std::string&, uint32_t&, std::string&, std::string&, std::vector<std::string>&, std::string&, std::string&, std::string& >() > start;
	};

} // end namespace graphite

*/
#endif //GRAPHITE_VCFPARSER_H
