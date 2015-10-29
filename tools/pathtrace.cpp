#include "core/util/Params.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/reference/FastaWriter.h"
#include "plugins/pathtrace/graph/VariantGraph.h"
#include "plugins/pathtrace/graph/SNPNode.h"

void printNodesToFasta(const std::string& filePrefix, std::vector< std::vector< graphite::INode::SharedPtr > >& nodes, const std::string& region);
void printNodes(std::vector< graphite::INode::SharedPtr >& nodes, std::ostream& out, size_t id, const std::string& region);

int main(int argc, char** argv)
{
	graphite::Params params;
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

	auto fastaReferencePtr = std::make_shared< graphite::FastaReference >(fastaPath, regionPtr);

	// load variants from vcf
	auto variantManagerPtr = std::make_shared< graphite::VCFManager >(vcfPaths, regionPtr);
	variantManagerPtr->asyncLoadVCFs(); // begin the process of loading the vcfs asynchronously
	variantManagerPtr->waitForVCFsToLoadAndProcess(); // wait for vcfs to load into memory
	variantManagerPtr->releaseResources(); // releases the vcf file memory, we no longer need the file resources
	auto variantListPtr = variantManagerPtr->getVariantsInRegion(regionPtr);

	auto variantGraphPtr = std::make_shared< graphite::vg::VariantGraph >(fastaReferencePtr, variantListPtr);
	variantGraphPtr->constructGraph();

	std::vector< std::string > paths;
	std::vector< std::vector< graphite::INode::SharedPtr > > nodes;
	variantGraphPtr->getAllPaths(paths, nodes);
	printNodesToFasta(filePrefix, nodes, regionPtr->getRegionString());

	return 0;
}

void printNodes(std::vector< graphite::INode::SharedPtr >& nodes, std::ostream& out, size_t id, const std::string& region)
{
	std::string pathString;
	std::string header = std::to_string(id) + " Region=" + region;
	std::string variantsString = "Variants=";
	uint32_t variantsCount = 0;
	for (auto node : nodes)
	{
		auto referenceNodePtr = std::dynamic_pointer_cast< graphite::vg::ReferenceNode >(node);
		auto variantNodePtr = std::dynamic_pointer_cast< graphite::vg::SNPNode >(node);
		if (variantNodePtr != nullptr)
		{
			std::string suffix = (variantsCount > 0) ? "," : "";
			variantsString += suffix + std::to_string(variantNodePtr->getPosition()) + ":" + variantNodePtr->getSequence() + ":" + variantNodePtr->getReferenceSequence();
			++variantsCount;
		}
		// header += std::string(node->getSequence(), node->getLength()) + "|";
		pathString += std::string(node->getSequence(), node->getLength());
	}
	if (variantsCount > 0)
	{
		header += " " + variantsString;
	}
	graphite::FastaWriter fastaWriter(header, pathString);
	fastaWriter.write(out);
}

void printNodesToFasta(const std::string& filePrefix, std::vector< std::vector< graphite::INode::SharedPtr > >& nodes, const std::string& region)
{
	size_t count = 1;
	for (auto nodeList : nodes)
	{
		if (filePrefix.size() > 0)
		{
			std::string fileName = filePrefix + "_" + std::to_string(count++) + ".fa";
			ofstream out;
			out.open(fileName, std::ios::out | std::ios::trunc);
			printNodes(nodeList, out, (count - 1), region);
			out.close();
		}
		else
		{
			printNodes(nodeList, std::cout, (count - 1), region);
			std::cout << std::endl;
		}
		++count;
	}
}
