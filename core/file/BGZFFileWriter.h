#ifndef GRAPHITE_BGZFFILEWRITER_H
#define GRAPHITE_BGZFFILEWRITER_H

#include "IFileWriter.h"

#include <iostream>
#include <stdint.h>

namespace graphite
{
	class BGZFFileWriter : public IFileWriter
	{
	public:
		typedef std::shared_ptr< BGZFFileWriter > SharedPtr;

		BGZFFileWriter(const std::string& filePath);
		~BGZFFileWriter();

		bool open() override;
		void close() override;
		bool write(const char* data, size_t dataLength) override;
		std::string getFilePath() override { return m_file_path; }

	private:
		int deflateBlock();
		void flushBlock();
		// 'packs' an unsigned integer into the specified buffer
		inline void packUnsignedInt(char* buffer, unsigned int value)
		{
			buffer[0] = (char)value;
			buffer[1] = (char)(value >> 8);
			buffer[2] = (char)(value >> 16);
			buffer[3] = (char)(value >> 24);
		}

		inline void packUnsignedShort(char* buffer, unsigned short value)
		{
			buffer[0] = (char)value;
			buffer[1] = (char)(value >> 8);
		}

		std::string m_file_path;
		uint32_t m_block_offset;
		size_t m_uncompressed_block_size;
		char* m_uncompressed_block;
		size_t m_compressed_block_size;
		char* m_compressed_block;
		uint64_t m_block_address;
		bool m_is_open;
		FILE* m_file;
		int m_window_size;

		// consts
		// zlib constants
		const int GZIP_ID1   = 31;
		const int GZIP_ID2   = 139;
		const int CM_DEFLATE = 8;
		const int FLG_FEXTRA = 4;
		const int OS_UNKNOWN = 255;
		const int BGZF_XLEN  = 6;
		const int BGZF_ID1   = 66;
		const int BGZF_ID2   = 67;
		const int BGZF_LEN   = 2;
		const int GZIP_WINDOW_BITS    = -15;
		const int Z_DEFAULT_MEM_LEVEL = 8;

		// BZGF constants
		const int BLOCK_HEADER_LENGTH = 18;
		const int BLOCK_FOOTER_LENGTH = 8;
		const int MAX_BLOCK_SIZE      = 65536;
		const int DEFAULT_BLOCK_SIZE  = 65536;
	};
}

#endif //GRAPHITE_BGZFFILEWRITER_H
