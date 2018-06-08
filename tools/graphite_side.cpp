#include <Python.h>
#include <string>
#include <iostream>

int main(int argc, char** argv)
{
	// get all param properties
	graphite::Params params;
	params.parseGSSW(argc, argv);
	if (params.showHelp() || !params.validateRequired())
	{
		params.printHelp();
		return 0;
	}
	auto bamPaths = params.getBAMPaths();
	auto fastaPath = params.getFastaPath();
	auto vcfPaths = params.getInVCFPaths();
	auto outputDirectory = params.getOutputDirectory();
	auto paramRegionPtr = params.getRegion();
	auto swPercent = params.getPercent();
	auto threadCount = params.getThreadCount();
	auto matchValue = params.getMatchValue();
	auto misMatchValue = params.getMisMatchValue();
	auto gapOpenValue = params.getGapOpenValue();
	auto gapExtensionValue = params.getGapExtensionValue();
	auto includeDuplicates = params.getIncludeDuplicates();
	auto graphSize = params.getGraphSize();
	auto outputVisualizationFiles = params.outputVisualizationFiles();

	// create reference reader
	// create VCF readers
	// create VCF writers
	// create bam readers
	// for each region:
	// create graph processor
	// call process on processor
	// close reference reader
	// close VCF writers
	// close bam readers
	// close VCF readers
	// exit

	/*
	// Py_SetProgramName(std::string("PyTest"));
	setenv("PYTHONPATH","/uufs/chpc.utah.edu/common/home/u0105808/Projects/graphite_side_graph/scripts/",1);
	PyObject *pName, *pModule, *pDict, *pFunc, *pValue, *presult, *pArgs;
	Py_Initialize();
	pName = PyString_FromString((char*)"parseVCF");
	pModule = PyImport_Import(pName);
	pFunc = PyObject_GetAttrString(pModule, "someFunction");
	if (PyCallable_Check(pFunc))
	{
		// auto pArg = PyString_FromString("jimmer");
		auto pArg=Py_BuildValue("(z)",(char*)"jimmer");
		pArgs = PyTuple_New(1);
		PyTuple_SetItem(pArgs, 0, pArg);
		presult = PyObject_CallObject(pFunc, pArg);
	}
	else
	{
		PyErr_Print();
	}
	std::string result = PyString_AsString(presult);
	std::cout << result << std::endl;
	// Py_XDECREF(repr);
	// Py_XDECREF(str);
	Py_Finalize();
	*/
	return 0;
}
