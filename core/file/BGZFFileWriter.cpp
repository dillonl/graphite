/*
 * Basically everything here has been lifted directly from Bamtools. Specifically this file:
 * http://mendel.stanford.edu/sidowlab/downloads/hot/UniPeak_1.0/misc/bamtools/BGZF.cpp
 *
 * Thanks Derek!
 */

#include "BGZFFileWriter.h"

#include "zlib.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

namespace graphite
{

	BGZFFileWriter::BGZFFileWriter(const std::string& filePath) :
		IFileWriter(FileType::BGZF),
		m_file_path(filePath),
		m_is_open(false),
		m_file(nullptr),
		m_uncompressed_block(nullptr),
		m_uncompressed_block_size(65536),
		m_compressed_block(nullptr),
		m_compressed_block_size(65536),
		m_block_offset(0),
		m_block_address(0)
	{
		m_uncompressed_block = new char[m_uncompressed_block_size];
		m_compressed_block = new char[m_compressed_block_size];
	}

	BGZFFileWriter::~BGZFFileWriter()
	{
		if(m_compressed_block) delete[] m_compressed_block;
		if(m_uncompressed_block) delete[] m_uncompressed_block;
	}

	bool BGZFFileWriter::open()
	{
		m_file = fopen(m_file_path.c_str(), "wb");
		if ( !m_file )
		{
			fprintf(stderr, "BGZF ERROR: unable to open file %s\n", m_file_path.c_str() );
			return false;
		}

		// set flags, return success
		m_is_open = true;
		return true;
	}

	void BGZFFileWriter::close()
	{
		// skip if file not open, otherwise set flag
		if ( !m_is_open) return;

		flushBlock();
		int blockLength = deflateBlock();
		fwrite(m_compressed_block, 1, blockLength, m_file);

		// flush and close
		fflush(m_file);
		fclose(m_file);
		m_is_open = false;
	}

	bool BGZFFileWriter::write(const char* data, size_t dataLength)
	{
		if (!m_is_open)
		{
			return false;
		}

		// initialize
		unsigned int numBytesWritten = 0;
		const char* input = data;
		unsigned int blockLength = m_uncompressed_block_size;

		// copy the data to the buffer
		while (numBytesWritten < dataLength)
		{

			unsigned int copyLength = std::min((size_t)(blockLength - m_block_offset), (size_t)(dataLength - numBytesWritten));
			char* buffer = m_uncompressed_block;
			memcpy(buffer + m_block_offset, input, copyLength);

			m_block_offset += copyLength;
			input += copyLength;
			numBytesWritten += copyLength;

			if ( m_block_offset == blockLength )
			{
				flushBlock();
			}
		}

		return numBytesWritten;
	}

	// flushes the data in the BGZF block
	void BGZFFileWriter::flushBlock(void)
	{

		// flush all of the remaining blocks
		while ( m_block_offset > 0 ) {

			// compress the data block
			int blockLength = deflateBlock();

			// flush the data to our output stream
			int numBytesWritten = fwrite(m_compressed_block, 1, blockLength, m_file);

			if (numBytesWritten != blockLength)
			{
				fprintf(stderr, "BGZF ERROR: expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
				exit(1);
			}

			m_block_address += blockLength;
		}
	}

	int BGZFFileWriter::deflateBlock()
	{
		// initialize the gzip header
		char* buffer = m_compressed_block;
		memset(buffer, 0, 18);
		buffer[0]  = GZIP_ID1;
		buffer[1]  = (char)GZIP_ID2;
		buffer[2]  = CM_DEFLATE;
		buffer[3]  = FLG_FEXTRA;
		buffer[9]  = (char)OS_UNKNOWN;
		buffer[10] = BGZF_XLEN;
		buffer[12] = BGZF_ID1;
		buffer[13] = BGZF_ID2;
		buffer[14] = BGZF_LEN;

		// set compression level
		const int compressionLevel = Z_DEFAULT_COMPRESSION;

		// loop to retry for blocks that do not compress enough
		int inputLength = m_block_offset;
		int compressedLength = 0;
		unsigned int bufferSize = m_compressed_block_size;

		while ( true )
		{

			// initialize zstream values
			z_stream zs;
			zs.zalloc    = NULL;
			zs.zfree     = NULL;
			zs.next_in   = (Bytef*)m_uncompressed_block_size;
			zs.avail_in  = inputLength;
			zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
			zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

			// initialize the zlib compression algorithm
			if (deflateInit2(&zs, compressionLevel, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK)
			{
				fprintf(stderr, "BGZF ERROR: zlib deflate initialization failed.\n");
				exit(1);
			}

			// compress the data
			int status = deflate(&zs, Z_FINISH);
			if ( status != Z_STREAM_END )
			{
				deflateEnd(&zs);

				// reduce the input length and try again
				if (status == Z_OK)
				{
					inputLength -= 1024;
					if(inputLength < 0)
					{
						fprintf(stderr, "BGZF ERROR: input reduction failed.\n");
						exit(1);
					}
					continue;
				}

				fprintf(stderr, "BGZF ERROR: zlib::deflateEnd() failed.\n");
				exit(1);
			}

			// finalize the compression routine
			if (deflateEnd(&zs) != Z_OK)
			{
				fprintf(stderr, "BGZF ERROR: zlib::deflateEnd() failed.\n");
				exit(1);
			}

			compressedLength = zs.total_out;
			compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
			if (compressedLength > MAX_BLOCK_SIZE)
			{
				fprintf(stderr, "BGZF ERROR: deflate overflow.\n");
				exit(1);
			}

			break;

			// store the compressed length
			packUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

			// store the CRC32 checksum
			unsigned int crc = crc32(0, NULL, 0);
			crc = crc32(crc, (Bytef*)m_uncompressed_block, inputLength);
			packUnsignedInt(&buffer[compressedLength - 8], crc);
			packUnsignedInt(&buffer[compressedLength - 4], inputLength);

			// ensure that we have less than a block of data left
			int remaining = m_block_offset - inputLength;
			if ( remaining > 0 )
			{
				if (remaining > inputLength)
				{
					fprintf(stderr, "BGZF ERROR: after deflate, remainder too large.\n");
					exit(1);
				}
				memcpy(m_uncompressed_block, m_uncompressed_block + inputLength, remaining);
			}

			m_block_offset = remaining;
			return compressedLength;
		}
	}
}
