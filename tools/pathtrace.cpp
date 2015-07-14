#include "core/util/Params.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/reference/FastaWriter.h"
#include "plugins/vg/graph/VariantGraph.h"

void printNodesToFasta(const std::string& filePrefix, std::vector< std::vector< gwiz::INode::SharedPtr > >& nodes);
void printNodes(std::vector< gwiz::INode::SharedPtr >& nodes, std::ostream& out);

int main(int argc, char** argv)
{
	gwiz::Params params;
	params.parsePathTrace(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		exit(0);
	}
	auto fastaPath = params.getFastaPath();
	auto vcfPaths = params.getInVCFPaths();
	auto regionPtr = params.getRegion();
	auto filePrefix = params.getFilePrefix();

	auto fastaReferencePtr = std::make_shared< gwiz::FastaReference >(fastaPath, regionPtr);

	// load variants from vcf
	auto variantManagerPtr = std::make_shared< gwiz::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	std::cout << "variants loaded" << std::endl;
	auto variantListPtr = variantManagerPtr->getVariantsInRegion(regionPtr);
	std::cout << "variants in region" << std::endl;

	auto variantGraphPtr = std::make_shared< gwiz::vg::VariantGraph >(fastaReferencePtr, variantListPtr);
	std::cout << "starting graph construction" << std::endl;
	variantGraphPtr->constructGraph();

	std::vector< std::string > paths;
	std::vector< std::vector< gwiz::INode::SharedPtr > > nodes;
	variantGraphPtr->getAllPaths(paths, nodes);
	printNodesToFasta(filePrefix, nodes);

	return 0;
}

void printNodes(std::vector< gwiz::INode::SharedPtr >& nodes, std::ostream& out)
{
	std::string pathString;
	for (auto node : nodes)
	{
		pathString += std::string(node->getSequence(), node->getLength());
	}
	std::string header = "test";
	gwiz::FastaWriter fastaWriter(header, pathString);
	fastaWriter.write(out);
}

void printNodesToFasta(const std::string& filePrefix, std::vector< std::vector< gwiz::INode::SharedPtr > >& nodes)
{
	size_t count = 1;
	for (auto nodeList : nodes)
	{
		if (filePrefix.size() > 0)
		{
			std::string fileName = filePrefix + "_" + std::to_string(count++);
			ofstream out;
			out.open(fileName, std::ios::out | std::ios::trunc);
			printNodes(nodeList, out);
			out.close();
		}
		else
		{
			printNodes(nodeList, std::cout);
		}
	}
}
