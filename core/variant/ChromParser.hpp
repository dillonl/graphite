#ifndef GWIZ_CHROMPARSER_H
#define GWIZ_CHROMPARSER_H

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
	struct ChromParser : boost::spirit::qi::grammar<Iterator, position() >
	{
		ChromParser() : ChromParser::base_type(start)
			{
				using namespace boost::spirit;

				tab = qi::lit('\t');
				chrom %= +qi::graph;
				pos %= qi::uint_;
				start %= chrom > tab > pos > tab;

				// names for error output
				tab.name("TAB");
				pos.name("POS");

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


		boost::spirit::qi::rule<Iterator, void()> tab;
		boost::spirit::qi::rule<Iterator, void()> chrom;
		boost::spirit::qi::rule<Iterator, position()> pos;
		boost::spirit::qi::rule<Iterator, position() > start;
	};

} // end namespace gwiz

#endif //GWIZ_CHROMPARSER_H
