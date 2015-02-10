#include "gtest/gtest.h"

#include "ssw/ssw_cpp.h"

namespace
{

	TEST(SSWTests, SSW_ALIGN)
	{
		const std::string ref   = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAACAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAACAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAACAGCCTTTCTGACCCGGAAATCAAAATAGGCAC";
		const std::string query = "ATTTTTTCGCCATTTCCATTTTTAAAAAAAAAAAAACCCCCCCCGGGGGGTTTTTTTTTTTTTTTTTTTTCACCCCCCCCCACACACACACACACACACACAAACACACACACACACACACAATTTTTTCGCCATTTCCATTTTT";
		//const string ref   = "CCGTTTATCGCA";
		//const string query = "CCTTTTATCGCA";

		for (auto i = 0; i < 1000000; ++i)
		{
			// Declares a default Aligner
			StripedSmithWaterman::Aligner aligner;
			// Declares a default filter
			StripedSmithWaterman::Filter filter;
			// Declares an alignment that stores the result
			StripedSmithWaterman::Alignment alignment;
			// Aligns the query to the ref
			aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);
		}

		// PrintAlignment(alignment);
	}

	TEST(SSWTests, SSW_ALIGN_ALL_MATCH)
	{
		const std::string ref   = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
		const std::string query = "CAGCCTTTCTGACCCGGAAATCAAAATAGGCACAACAAA";
		//const string ref   = "CCGTTTATCGCA";
		//const string query = "CCTTTTATCGCA";

		// Declares a default Aligner
		StripedSmithWaterman::Aligner aligner;
		// Declares a default filter
		StripedSmithWaterman::Filter filter;
		// Declares an alignment that stores the result
		StripedSmithWaterman::Alignment alignment;
		// Aligns the query to the ref
		aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment);

		PrintAlignment(alignment);
	}

}
