#ifndef GWIZ_VCFPARSER_H
#define GWIZ_VCFPARSER_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <utility>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>


namespace gwiz
{

	template<typename Iterator>
	struct VariantParser : boost::spirit::qi::grammar<Iterator, boost::fusion::vector< std::string&, uint32_t&, std::string&, std::vector<std::string>&, std::vector<std::string>&, std::string&, std::string&, std::unordered_map< std::string, std::string >& >() >
	{
		VariantParser() : VariantParser::base_type(start)
			{
				using namespace boost::spirit;

				std::string refAltMatchString = "actgnACTGN";
				auto refAltMatch = (+qi::char_(refAltMatchString) >> *(',' >> +qi::char_(refAltMatchString))) || +qi::graph;
				tab = qi::lit('\t');
				chrom %= +qi::graph;
				pos %= qi::uint_;
				id %= +qi::graph;
				ref %= refAltMatch;
				alt %= refAltMatch;
				qual = +qi::graph;
				filter = +qi::graph;

				info = pair >> *((qi::lit(';') | '&') >> pair);
				pair = key >> -('=' >> value);
				key = qi::char_("a-zA-Z_0-9\\./(),") >> *qi::char_("a-zA-Z_0-9\\./(),");
				value = +qi::char_("a-zA-Z_0-9\\./(),");
				start %= chrom > tab > pos > tab > id > tab > ref > tab > alt > tab > qual > tab > filter > tab > info;

				// names for error output
				tab.name("TAB");
				chrom.name("CHROM");
				pos.name("POS");
				id.name("ID");
				ref.name("REF");
				alt.name("ALT");
				qual.name("QUAL");
				filter.name("FILTER");
				info.name("INFO");

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
		boost::spirit::qi::rule<Iterator, std::string()> chrom;
		boost::spirit::qi::rule<Iterator, uint32_t()> pos;
		boost::spirit::qi::rule<Iterator, std::string()> id;
		boost::spirit::qi::rule<Iterator, std::vector<std::string>() > ref;
		boost::spirit::qi::rule<Iterator, std::vector<std::string>() > alt;
		boost::spirit::qi::rule<Iterator, std::string()> qual;
		boost::spirit::qi::rule<Iterator, std::string()> filter;
		boost::spirit::qi::rule<Iterator, std::unordered_map< std::string, std::string >() > info;
		boost::spirit::qi::rule<Iterator, std::pair< std::string, std::string >() > pair;
		boost::spirit::qi::rule<Iterator, std::string() > key, value;

		boost::spirit::qi::rule<Iterator, boost::fusion::vector< std::string&, uint32_t&, std::string&, std::vector<std::string>&, std::vector<std::string>&, std::string&, std::string&, std::unordered_map< std::string, std::string >& >() > start;
			// boost::spirit::qi::rule<Iterator, boost::fusion::vector< std::string&, uint32_t&, std::string&, std::vector<std::string>&, std::vector<std::string>&, std::string&, std::string& >() > start;
	};

} // end namespace gwiz

#endif //GWIZ_VCFPARSER_H
