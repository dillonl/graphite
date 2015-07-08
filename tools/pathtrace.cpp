#include "core/util/Params.h"
#include "core/variant/VCFManager.h"
#include "core/reference/FastaReference.h"
#include "core/reference/FastaWriter.h"
#include "plugins/vg/graph/VariantGraph.h"

void printNodesToFasta(std::vector< std::vector< gwiz::INode::SharedPtr > >& nodes, std::ostream out);

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
	auto outputVCFPath = params.getOutVCFPath();
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

	return 0;
}

void printNodesToFasta(std::vector< std::vector< gwiz::INode::SharedPtr > >& nodes, std::ostream out)
{
	/*
	  size_t count = 1;
	  for (auto nodeList : nodes)
	  {
	  std::ostream out;
	  if (filePrefix.size() > 0)
	  {
	  ofstream outfileStream;
	  out = outFileStream;
	  out.open(fileName = "" + std::to_string(count++), std::ios::out | std::ios::trunc);
	  }
	  else
	  {
	  out = std::cout;
	  }
	  std::string pathString;
	  for (auto node : nodePaths)
	  {
	  pathString += std::string(node->getSequence(), node->getLength());
	  }
	  std::string header = "test";
	  gwiz::FastaWriter fastaWriter(header, pathString);
	  fastaWriter.write(out);
	  }
	*/
}
