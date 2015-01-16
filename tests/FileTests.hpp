#ifndef GWIZ_TESTS_FILETESTS_HPP
#define GWIZ_TESTS_FILETESTS_HPP

#include "utils/file/ASCIIFileReader.h"
#include "TestConfig.h"

TEST(ASCIIFileReaderTests, OpenValidFileTest)
{
	// auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader::SharedPtr>(TEST_LINE_NUMBERS_FILE);
	// auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>("");
	bool success = true;
	try
	{
		std::string path = TEST_LINE_NUMBERS_FILE;
		auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
		asciiReaderPtr->Open();
	}
	catch (std::ios_base::failure& exc)
	{
		success = false;
	}
	EXPECT_TRUE(success);
}

TEST(ASCIIFileReaderTests, OpenInvalidFileTest)
{
	bool success = true;
	try
	{
		std::string path = TEST_INVALID_FILE;
		auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
		asciiReaderPtr->Open();
	}
	catch (std::ios_base::failure& exc)
	{
		success = false;
	}
	EXPECT_FALSE(success);
}

TEST(ASCIIFileReaderTests, getNextLineFileTest)
{
	bool success = true;
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
	asciiReaderPtr->Open();
	const char* line;
	uint32_t count = 1;
	while ((line = asciiReaderPtr->getNextLine()) != NULL)
	{
		ASSERT_EQ(atoi(line), count);
		++count;
	}
	EXPECT_TRUE(success);
}

TEST(ASCIIFileReaderTests, getNextLineEndIsNULLFileTest)
{
	std::string path = TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
	asciiReaderPtr->Open();
	for (uint32_t i = 1; i <= 100; ++i)
	{
		asciiReaderPtr->getNextLine();
	}
	const char* line = asciiReaderPtr->getNextLine();
	bool success = line == NULL;
	EXPECT_TRUE(success);
}

TEST(ASCIIFileReaderTests, getNextLineUnopenedFileTest)
{
	std::string path =  TEST_LINE_NUMBERS_FILE;
	auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
	bool success = (asciiReaderPtr->getNextLine() == NULL);
	EXPECT_TRUE(success);
}

/*

TEST(ASCIIFileReaderTests, CountLinesFileTest)
{
	bool success = true;
	try
	{
		std::string path = "/home/dlee/data/var.vcf";
		auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
		uintmax_t lineCount = asciiReaderPtr->CountLines(path);
		std::cout << "--lines: " << lineCount << std::endl;
	}
	catch (std::ios_base::failure& exc)
	{
		success = false;
	}
	EXPECT_TRUE(success);
}
/*
TEST(ASCIIFileReaderTests, CountLines2FileTest)
{
	bool success = true;
	try
	{
		std::string path = "/home/dlee/data/var.vcf";
		auto asciiReaderPtr = std::make_shared<gwiz::ASCIIFileReader>(path);
		uintmax_t lineCount = asciiReaderPtr->CountLines2(path);
	}
	catch (std::ios_base::failure& exc)
	{
		success = false;
	}
	EXPECT_TRUE(success);
}
*/
#endif //GWIZ_TESTS_FILETESTS_HPP
