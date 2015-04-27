#ifndef GWIZ_TESTS_FILETESTS_HPP
#define GWIZ_TESTS_FILETESTS_HPP

#include "core/utils/file/ASCIIFileReader.h"
#include "core/utils/file/ASCIIGZFileReader.h"
#include "TestConfig.h"

TEST(ASCIIGZFileReaderTests, OpenValidFileTest)
{
	bool success = true;
	try
	{
		std::string path = TEST_LINE_NUMBERS_GZ_FILE;
		auto asciiGZReaderPtr = std::make_shared<gwiz::ASCIIGZFileReader>(path);
		asciiGZReaderPtr->Open();
	}
	catch (std::ios_base::failure& exc)
	{
		success = false;
	}
	EXPECT_TRUE(success);
}

TEST(ASCIIGZFileReaderTests, getNextLineFileTest)
{
	bool success = false;
	std::string path = TEST_LINE_NUMBERS_GZ_FILE;
	auto asciiReaderPtr = std::make_shared<gwiz::ASCIIGZFileReader>(path);
	asciiReaderPtr->Open();
	const char* line;
	uint32_t count = 1;
	asciiReaderPtr->getNextLine();

	/*
	while ((line = asciiReaderPtr->getNextLine()) != nullptr)
	{
		std::string countString = std::to_string(count);
		ASSERT_STREQ(line, countString.c_str());
		++count;
		success = true;
	}
	EXPECT_TRUE(success);
	*/
}

/*
TEST(ASCIIFileReaderTests, OpenValidFileTest)
{
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
