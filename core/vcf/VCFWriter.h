

#ifndef GRAPHITE_VCFWRITER_H
#define GRAPHITE_VCFWRITER_H

#include "core/util/Noncopyable.hpp"
#include "core/sample/Sample.h"

#include <memory>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

namespace graphite
{
	class VCFWriter : private Noncopyable
	{
	public:
		typedef std::shared_ptr< VCFWriter > SharedPtr;
		VCFWriter(const std::string& filename, std::vector< graphite::Sample::SharedPtr >& bamSamplePtrs, const std::string& outputDirectory);
		~VCFWriter();

		void writeLine(const std::string& line);
        void writeHeader(const std::vector< std::string >& headerLines);
		Sample::SharedPtr getSamplePtr(const std::string& sampleName);
		std::vector< std::string > getSampleNames();
		std::vector< std::string > getColumnNames();
		bool isSampleNameInOriginalVCF(const std::string& sampleName);
		bool isSampleNameInBam(const std::string& sampleName);
		void setBlankFormatString(const std::string& blankFormatString);
		std::shared_ptr< std::string > getBlankFormatStringPtr();

	private:
		std::vector< Sample::SharedPtr > m_bam_sample_ptrs;
		std::unordered_map< std::string, bool > m_sample_name_in_vcf;
		std::vector< std::string > m_vcf_column_names;
		std::vector< std::string > m_sample_names;
		std::ofstream m_out_file;
		std::shared_ptr< std::string > m_black_format_string;
		std::unordered_map< std::string, Sample::SharedPtr > m_bam_sample_ptrs_map;
		std::vector< std::tuple< std::string, std::string > > m_format = {std::make_tuple("ID=DP_NFP", "##FORMAT=<ID=DP_NFP,Number=1,Type=Integer,Description=\"Read count at 95 percent Smith Waterman score or above\">"),
																		  std::make_tuple("ID=DP_NP", "##FORMAT=<ID=DP_NP,Number=1,Type=Integer,Description=\"Read count between 90 and 94 percent Smith Waterman score\">"),
																		  std::make_tuple("ID=DP_EP", "##FORMAT=<ID=DP_EP,Number=1,Type=Integer,Description=\"Read count between 80 and 89 percent Smith Waterman score\">"),
																		  std::make_tuple("ID=DP_SP", "##FORMAT=<ID=DP_SP,Number=1,Type=Integer,Description=\"Read count between 70 and 79 percent Smith Waterman score\">"),
																		  std::make_tuple("ID=DP_LP", "##FORMAT=<ID=DP_LP,Number=1,Type=Integer,Description=\"Read count at 69 percent or less Smith Waterman score\">"),
																		  std::make_tuple("ID=DP_AP", "##FORMAT=<ID=DP_AP,Number=1,Type=Integer,Description=\"Read count for mappings which map equally well into (or out of) reference and variant. Not resolvable but valid mapping.\">"),
																		  std::make_tuple("ID=DP4_NFP", "##FORMAT=<ID=DP4_NFP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling at 95 percent Smith Waterman score or above.\">"),
																		  std::make_tuple("ID=DP4_NP", "##FORMAT=<ID=DP4_NP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 90 and 94 percent Smith Waterman score.\">"),
																		  std::make_tuple("ID=DP4_EP", "##FORMAT=<ID=DP4_EP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 80 and 89 percent Smith Waterman score.\">"),
																		  std::make_tuple("ID=DP4_SP", "##FORMAT=<ID=DP4_SP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling between 70 and 79 percent Smith Waterman score.\">"),
																		  std::make_tuple("ID=DP4_LP", "##FORMAT=<ID=DP4_LP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling at 69 percent or less Smith Waterman score.\">"),
		                                                                  std::make_tuple("ID=DP4_AP", "##FORMAT=<ID=DP4_AP,Number=.,Type=Integer,Description=\"Number of 1) forward ref alleles, 2) reverse ref, 3) forward non-ref, 4) reverse non-ref alleles, used in variant calling for reads which map equally well into (or out of) reference and the alt. Not resolvable but a valid mapping.\">")};
	};
}

#endif //GRAPHITE_VCFWRITER_H
