#ifndef GWIZ_REGIONPARSER_H
#define GWIZ_REGIONPARSER_H

#include "core/util/Types.h"

#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>

namespace gwiz
{

	template<typename Iterator>
	struct RegionParser : boost::spirit::qi::grammar<Iterator, boost::fusion::vector< std::string&, position&, position& >() >
	{
		RegionParser() : RegionParser::base_type(start)
			{
				using namespace boost::spirit;

				hyphen = qi::lit('-');
				chrom = *~qi::char_(':');
				start_position %= qi::uint_;
				end_position %= qi::uint_;
				start %= chrom >> qi::lit(':') >> start_position >> hyphen >> end_position;

				// names for error output
				hyphen.name("HYPHEN");
				start_position.name("START_POSITION");
				end_position.name("END_POSITION");

				qi::on_error<qi::fail>
					(
						start
						, std::cout
						<< boost::phoenix::val("Error! Column ")
						<< _4                               // what failed?
						<< boost::phoenix::val(" is incorrectly formatted: \"")
						<< _1
						<< std::endl
						);

			}


		boost::spirit::qi::rule<Iterator, void()> hyphen;
		boost::spirit::qi::rule<Iterator, std::string()> chrom;
		boost::spirit::qi::rule<Iterator, position()> start_position;
		boost::spirit::qi::rule<Iterator, position()> end_position;
		boost::spirit::qi::rule<Iterator, boost::fusion::vector< std::string&, position&, position& >() > start;
	};

} // end namespace gwiz

#endif //GWIZ_REGIONPARSER_H
