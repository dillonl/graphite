#ifndef GRAPHITE_INFOFIELDPARSER_H
#define GRAPHITE_INFOFIELDPARSER_H

/*
#include "core/util/Types.h"

#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

namespace graphite
{

	namespace qi = boost::spirit::qi;

	template <typename Iterator>
	struct InfoFieldParser
		: qi::grammar<Iterator, std::map<std::string, std::string>()>
	{
		InfoFieldParser()
			: InfoFieldParser::base_type(query)
			{
				query =  pair >> *((qi::lit(';')) >> pair);
				pair  =  key >> -('=' >> value);
				key   =  qi::char_("a-zA-Z_") >> *qi::char_("a-zA-Z_0-9");
				value = +(qi::char_(",a-zA-Z_0-9-"));
			}
		qi::rule<Iterator, std::map<std::string, std::string>()> query;
		qi::rule<Iterator, std::pair<std::string, std::string>()> pair;
		qi::rule<Iterator, std::string()> key, value;
	};


} // end namespace graphite

*/
#endif //GRAPHITE_INFOFIELDPARSER_H
