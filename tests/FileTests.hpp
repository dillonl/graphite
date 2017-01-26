#ifndef GRAPHITE_TESTS_FILETESTS_HPP
#define GRAPHITE_TESTS_FILETESTS_HPP

#include "core/file/ASCIIFileReader.h"
#include "core/file/ASCIIGZFileReader.h"
#include "TestConfig.h"

namespace fileTests
{
	void testOpenValidFileTest(graphite::IFile::SharedPtr asciiFilePtr)
	{
		ASSERT_NO_THROW(asciiFilePtr->Open());
	}

	  // for some reason this causes the Travis CI to crash so I'll comment
	  // it out for now.
	void testOpenInvalidFileTest(graphite::IFile::SharedPtr asciiFilePtr)
	{
		ASSERT_THROW(asciiFilePtr->Open(), std::ios_base::failure);
	}

	void testReadNextLineUnopenedFile(graphite::IFile::SharedPtr asciiFilePtr)
	{
		std::string line;
		ASSERT_FALSE(asciiFilePtr->getNextLine(line));
	}

	void testReadNextLine(graphite::IFile::SharedPtr asciiFilePtr)
	{
		asciiFilePtr->Open();

		uint32_t count = 0;
		std::string line;
		while (asciiFilePtr->getNextLine(line))
		{
			++count;
			std::string countString = std::to_string(count);
			ASSERT_STREQ(line.c_str(), countString.c_str());
		}
		ASSERT_EQ(count, 100); // 100 is the last line of the text file
	}

	void testReadNextLineEOFReturnFalse(graphite::IFile::SharedPtr asciiFilePtr)
	{
		asciiFilePtr->Open();
		std::string line;
		for (uint32_t i = 1; i <= 100; ++i)
		{
			asciiFilePtr->getNextLine(line);
		}
		EXPECT_FALSE(asciiFilePtr->getNextLine(line));
	}
}

TEST(ASCIIGZFileReaderTests, OpenValidFileTest)
{
	std::string path = TEST_LINE_NUMBERS_GZ_FILE;
	std::cout << path << std::endl;
	auto gzAsciiReaderPtr = std::make_shared<graphite::ASCIIGZFileReader>(path);
	fileTests::testOpenValidFileTest(gzAsciiReaderPtr);
}

TEST(ASCIIFileReaderTests, OpenValidFileTest)
{
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<graphite::ASCIIFileReader>(path);
	fileTests::testOpenValidFileTest(asciiReaderPtr);
}

TEST(ASCIIGZFileReaderTests, ReadNextLineGZTest)
{
	std::string path = TEST_LINE_NUMBERS_GZ_FILE;
	auto gzAsciiReaderPtr = std::make_shared<graphite::ASCIIGZFileReader>(path);
	fileTests::testReadNextLine(gzAsciiReaderPtr);
}

TEST(ASCIIFileReaderTests, ReadNextLineTest)
{
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<graphite::ASCIIFileReader>(path);
	fileTests::testReadNextLine(asciiReaderPtr);
}

TEST(ASCIIGZFileReaderTests, ReadNextLineGZEOFReturnFalseTest)
{
	std::string path = TEST_LINE_NUMBERS_GZ_FILE;
	auto gzAsciiReaderPtr = std::make_shared<graphite::ASCIIGZFileReader>(path);
	fileTests::testReadNextLineEOFReturnFalse(gzAsciiReaderPtr);
}

TEST(ASCIIFileReaderTests, ReadNextLineEOFReturnFalseTest)
{
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<graphite::ASCIIFileReader>(path);
	fileTests::testReadNextLineEOFReturnFalse(asciiReaderPtr);
}

/*


TEST(ASCIIGZFileReaderTests, OpenInvalidFileTest)
{
std::string path = TEST_INVALID_FILE;
auto gzAsciiReaderPtr = std::make_shared<graphite::ASCIIGZFileReader>(path);
fileTests::testOpenInvalidFileTest(gzAsciiReaderPtr);
}

TEST(ASCIIFileReaderTests, OpenInvalidFileTest)
{
std::string path = TEST_INVALID_FILE;
auto asciiReaderPtr = std::make_shared<graphite::ASCIIFileReader>(path);
fileTests::testOpenInvalidFileTest(asciiReaderPtr);
}

TEST(ASCIIGZFileReaderTests, ReadNextLineGZUnopened)
{
	std::string path = TEST_LINE_NUMBERS_GZ_FILE;
	auto gzAsciiReaderPtr = std::make_shared<graphite::ASCIIGZFileReader>(path);
	fileTests::testReadNextLineUnopenedFile(gzAsciiReaderPtr);
}

TEST(ASCIIFileReaderTests, ReadNextLineUnopened)
{
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<graphite::ASCIIFileReader>(path);
	fileTests::testReadNextLineUnopenedFile(asciiReaderPtr);
}
*/

#endif //GRAPHITE_TESTS_FILETESTS_HPP
